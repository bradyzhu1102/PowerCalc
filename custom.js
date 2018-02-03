// added package from original
// Normal.js
// OnePCIGUI.js
// OneTCIGUI.js
// UniFunction.js
// Solve.js
// pc.js
// T.js
// Beta.js
// MoreMath.js
// Binomial.js
// PowerCalculator.js
// OnePGUI.js
// OneTGUI.js
// Chi2.js
// Tukey.js
// Poisson.js
// Pilot.js
// Rsquare.js
// RsquareGUI.js
// SimpleChi2GUI.js
// SimplePoissonGUI.js
// TwoTGUI.js


//begin of Normal.js
var NormalAux = {
  xMax : 9.0E+300,
  xMin : -9.0E+300,
  closedMax : false,
  closedMin : false,
  maxIter : 100,
  maxSearch : 25,
  xeps : 1.0E-006,
  feps : 1.0E-010,
  verbose : false
}
NormalAux.of = function(param){
   return Normal.cdf(param);
}
var Normal = {};
function NaNs(){ return 0.0/0.0;}
Normal.quantile = function(paramDouble){
    if ((paramDouble <= 0.0) || (paramDouble >= 1.0)) {
         return 0.0;
       }
    var d = 4.91 * (Math.pow(paramDouble, 0.14) - Math.pow(1.0 - paramDouble, 0.14));
    var x = Normal.search(NormalAux, paramDouble, d - 0.0025, 0.005);
    return x;
}
Normal.quantile3 = function(paramDouble1, paramDouble2, paramDouble3){
   var x = paramDouble2 + paramDouble3 * this.quantile(paramDouble1);
  return x;
}

Normal.search = function(paramUni,paramDouble1,paramDouble2, paramDouble3){
    var d1 = paramDouble2 + paramDouble3;
       var i = 0;

       if ((paramDouble2 + paramDouble1 + paramDouble3)==null ) {
         return NaNs();
       }
       var d2 = paramUni.of(paramDouble2) - paramDouble1;
       var d3 = paramUni.of(d1) - paramDouble1;
       if (d3==null) {
         return NaNs();
       }
       while (d2 * d3 > 0.0){
         i++;
         if (i > paramUni.maxSearch) {
           return NaNs();
         }
         paramDouble3 = -1.5 * (d1 - paramDouble2) * d3 / (d3 - d2);
         
         paramDouble2 = d1;
         d2 = d3;
         d1 += paramDouble3;
         if (d1 <= paramUni.xMin) {
           d1 = paramUni.closedMin ? paramUni.xMin : paramUni.xMin + 0.1 * (paramDouble2 - paramUni.xMin);
         }
         if (d1 >= paramUni.xMax) {
           d1 = paramUni.closedMax ? paramUni.xMax : paramUni.xMax - 0.1 * (paramUni.xMax - paramDouble2);
         }
         
         d3 = Normal.cdf(d1) - paramDouble1;
         if (d3==null) {
           return NaNs();
         }
       }
       var x = Normal.illinois(paramUni, paramDouble1, paramDouble2, d1, d2, d3);
       return x;
}

Normal.cdf = function(paramDouble){
       var d1 = 0.2316419;var d2 = 0.31938153;var d3 = -0.356563782;
       var d4 = 1.781477937;var d5 = -1.821255978;var d6 = 1.330274429;
       if (paramDouble > 0.0) {
         var x = 1.0 - Normal.cdf(-paramDouble);
         return x;
       }
       var d7;
       if (paramDouble > -5.0)
       {
         d7 = 1.0 / (1.0 - d1 * paramDouble);
         d7 = Normal.pdf(paramDouble) * d7 * (d2 + d7 * (d3 + d7 * (d4 + d7 * (d5 + d7 * d6))));
       }
       else
       {
         d7 = paramDouble * paramDouble;
         d7 = Normal.pdf(paramDouble) * (1.0 - (1.0 - 3.0 * (1.0 - 5.0 * (1.0 - 7.0 / d7) / d7) / d7) / d7) / -paramDouble;
       }
       return d7;
}

Normal.cdf3 = function(paramDouble1, paramDouble2, paramDouble3){
   var x = this.cdf((paramDouble1 - paramDouble2) / paramDouble3);
  return x;
}

Normal.pdf = function(paramDouble){
   var x = 0.39894228 * Math.exp(-0.5 * paramDouble * paramDouble);
  return x;
}

Normal.illinois = function(paramUniFunction, paramDouble1, paramDouble2, paramDouble3, paramDouble4, paramDouble5){
  var i = 0;
  
  if ((paramDouble4 + paramDouble5)==null) {
    return NaNs();
  }
  if (paramDouble4 * paramDouble5 > 0.0) {
    return NaNs();
  }
  while (i++ < paramUniFunction.maxIter)
  {
    if (Math.abs(paramDouble4 - paramDouble5) <= 1.0E-099) {
      return NaNs();
    }
    var d1 = paramDouble3 - paramDouble5 * (paramDouble3 - paramDouble2) / (paramDouble5 - paramDouble4);
    
    var d2 = paramUniFunction.of(d1) - paramDouble1;
    if (d2==null) {
      return NaNs();
    }
    if (d2 * paramDouble5 > 0.0){
      paramDouble4 /= 2.0;
    }else{
      paramDouble2 = paramDouble3;
      paramDouble4 = paramDouble5;
    }
    paramDouble3 = d1;
    paramDouble5 = d2;
    if ((Math.abs(paramDouble2 - paramDouble3) <= paramUniFunction.xeps) || (Math.abs(d2) <= paramUniFunction.feps)) {
      return d1;
    }
  }
  return NaNs();
}

Normal.power = function(paramDouble1, paramInt, paramDouble2, paramDouble3){
  if (paramDouble3 <= 0.0) {
    return NaNs();
  }
  if ((paramDouble2 <= 0.0) || (paramDouble2 >= 1.0)) {
    return NaNs();
  }
  if (paramInt < 0) {
    return this.power(-paramDouble1, 1, paramDouble2, paramDouble3);
  }
  if (paramInt == 0)
  {
    var d = this.quantile(paramDouble2 / 2.0);
    return this.cdf3(d, paramDouble1, paramDouble3) + 1.0 - this.cdf3(-d, paramDouble1, paramDouble3);
  }
  var d = -this.quantile(paramDouble2);
  var x = 1.0 - this.cdf3(d, paramDouble1, paramDouble3);
  return x;
}

Normal.power3 = function(paramDouble1, paramInt, paramDouble2){
   var x = this.power(paramDouble1, paramInt, paramDouble2, 1.0);
  return x;
}
Normal.power1 = function(paramDouble){
   var x = this.power(saveMu1, saveTail, paramDouble);
  return x;
}
//end of Normal.js

//begin of OnePCIGUI.js
var onePCIGUI = {
      Sigma:0.51,zCrit:1.96
      };

onePCIGUI.gui = function(){
   onePCIGUI.isFinite=1;
   onePCIGUI.N=1000.0;
   onePCIGUI.worstCase=1;
   onePCIGUI.pi=0.5;
   onePCIGUI.ME=0.09297;
   onePCIGUI.n = 100.0;
   onePCIGUI.conf = 0.95;
}

onePCIGUI.news = function(){
   onePCIGUI.Sigma=0.51;
   onePCIGUI.zCrit=1.96;
   onePCIGUI.isFinite=0;
   onePCIGUI.N=0.0;
   onePCIGUI.worstCase=0;
   onePCIGUI.pi=0.0;
   onePCIGUI.ME=0.0;
   onePCIGUI.n = 0.0;
   onePCIGUI.conf = 0.0;
}

onePCIGUI.click = function(){
   this.zCrit = Normal.quantile(1.0 - (1.0 - this.conf) / 2.0);
   this.worstCase = (Math.abs(this.pi - 0.5) < 1.0E-012 ? 1 : 0);
   this.Sigma =Math.sqrt(this.pi * (1.0 - this.pi));
   this.ME =Math.max(this.ME, 0.001 * this.Sigma);
   for (var i = 0; i < 3; i++) {
      this.n = (this.zCrit * this.Sigma / this.ME);
      this.n *= this.n;
      if (this.isFinite == 1) {
         this.n = (1.0 + this.n / (1.0 + this.n / this.N));
      }
      this.n = Math.max(2.0, this.n);
   }
}

onePCIGUI.worstCase_changed = function() {
   if (this.worstCase == 1) {
      this.pi = 0.5;
   }
   this.click();
}

onePCIGUI.n_changed = function() {
   this.N = Math.max(2.0, Math.round(this.N));
   this.n = Math.max(2.0, Math.round(this.n));
   this.zCrit = Normal.quantile(1.0 - (1.0 - this.conf) / 2.0);
   this.ME = (this.zCrit * this.Sigma / Math.sqrt(this.n - 1.0));
   if (this.isFinite == 1) {
      this.ME *= Math.sqrt(1.0 - this.n / this.N);
   }
}
function StringNaN(o){
   if(o!=null && o!="" && o!=undefined)
      return false;
   return true;
}
module.exports.cfap = function(id, key, value){
   if(onePCIGUI.n==undefined){
      onePCIGUI.news();
      onePCIGUI.gui();
      onePCIGUI.click();
   }
   var mm = 0;
   if(id!=null && id!=''){
      if(key!=null && key!=""){
         if(key=="isFinite"){
           if(StringNaN(value)){
              onePCIGUI.isFinite=0;
           }else{
              onePCIGUI.isFinite=1;
           }
           onePCIGUI.click();
        }
        if(key=="worstCase"){
           if(StringNaN(value)){
              onePCIGUI.worstCase=0;
           }else{
              onePCIGUI.worstCase=1;
           }
           onePCIGUI.worstCase_changed();
        }
        if(key=="N"){
           if(!StringNaN(value)){
              onePCIGUI.N=parseFloat(value);
              onePCIGUI.click();
           }
        }
        if(key=="n"){
           if(!StringNaN(value)){
              onePCIGUI.n=parseFloat(value);
              onePCIGUI.n_changed();
           }
        }
        if(key=="ME"){
           if(!StringNaN(value)){
              onePCIGUI.ME = parseFloat(value);
              onePCIGUI.click();
           }
        }
        if(key=="pi"){
           if(!StringNaN(value)){
              onePCIGUI.pi = parseFloat(value);
              onePCIGUI.click();
           }
        }
        if(key=="conf"){
           if(!StringNaN(value)){
              onePCIGUI.conf = parseFloat(value);
              onePCIGUI.click();
           }
        }
     }
  }else{
  onePCIGUI.news();
  onePCIGUI.gui();   
  }
  return onePCIGUI;
  }
//end of OnePCIGUI.js

//begin of OneTCIGUI.js
var OneTCIGUI = {}

OneTCIGUI.init = function(conf, nN, n, sigma, mE, tCrit, isFinite) {
  this.conf = conf;
  this.NN = nN;
  this.n = n;
  this.Sigma = sigma;
  this.ME = mE;
  this.tCrit = tCrit;
  OneTCIGUI.isFinite = isFinite;
}

OneTCIGUI.news = function() {
  this.conf = 0.0;
  this.NN = 0.0;
  this.n = 0.0;
  this.Sigma = 0.0;
  this.ME = 0.0;
  this.tCrit = 0.0;
  this.isFinite = 0;
}

OneTCIGUI.gui = function(){
  this.isFinite=1; 
  this.NN=1000.0;
  this.conf=0.95;
  this.Sigma=1.0;
  this.ME=0.2;
  this.n=25.0;    
  OneTCIGUI.n_changed();
}

OneTCIGUI.click = function(){
  this.tCrit = Normal.quantile(1.0 - (1.0 - this.conf) / 2.0);
  this.ME =Math.max(this.ME, 0.001 * this.Sigma);
  for (var i = 0; i < 3; i++){
    this.n = (this.tCrit * this.Sigma / this.ME);
    this.n *= this.n;
    if (this.isFinite == 1) {
      this.n /= (1.0 + this.n / this.NN);
    }
    this.n = Math.max(2.0, this.n);
    this.tCrit = T.quantile2(1.0 - (1.0 - this.conf) / 2.0, this.n - 1.0);
  }
}

OneTCIGUI.n_changed = function(){
  this.NN =Math.max(2.0, Math.round(this.NN));
  this.n = Math.max(2.0, Math.round(this.n));
  this.tCrit = T.quantile2(1.0 - (1.0 - this.conf) / 2.0, this.n - 1.0);
  this.ME = (this.tCrit * this.Sigma / Math.sqrt(this.n));
  if (this.isFinite == 1) {
    this.ME *= Math.sqrt(1.0 - this.n / this.NN);
  }
}

OneTCIGUI.isFinite_changed = function(){
  pc.setVisible("N", this.isFinite == 1);
  click();
}

module.exports.CImean = function(id, key, value){
  //console.log("CImean hit");
    if(isNullOrEmpty(OneTCIGUI.conf)){
      //console.log("CImean init hit");
     OneTCIGUI.news();
      OneTCIGUI.gui();
      OneTCIGUI.click();
    }
    //console.log(OneTCIGUI.isFinite);
    if(!isNullOrEmpty(id)){
          if(!isNullOrEmpty(key)){
            if(key=="isFinite"){
              if(isNullOrEmpty(value)){
                OneTCIGUI.isFinite=0;
              }else{
                OneTCIGUI.isFinite = parseInt(value);
              }
              OneTCIGUI.click();
            }
            if(key=="NN"){
              if(!isNullOrEmpty(value)){
                OneTCIGUI.NN=parseFloat(value);
              }
              OneTCIGUI.click();
            }
            if(key=="conf"){
               if(!isNullOrEmpty(value)){
                  OneTCIGUI.conf = parseFloat(value);
               }
               OneTCIGUI.click();
            }
            if(key=="sgn"){
               if(!isNullOrEmpty(value)){
                 OneTCIGUI.Sigma = parseFloat(value);
               }
               OneTCIGUI.click();
            }
            if(key=="me"){
               if(!isNullOrEmpty(value)){
                 OneTCIGUI.ME = parseFloat(value);
               }
               OneTCIGUI.click();
            }
            if(key=="n"){
              if(!isNullOrEmpty(value)){
                OneTCIGUI.n = parseFloat(value);
                 OneTCIGUI.n_changed();
              }
            }
          }
      }
      return OneTCIGUI;
}
//end of OneTCIGUI.js

//begin of UniFunction.js
function isNullOrEmpty(value){
  if(value!=null && value!=undefined && value!='')
    return false;
  return true;
}

function $G(){
  var Url=window.location.href;//如果想获取框架顶部的url可以用 top.window.location.href
  var u,g,StrBack='';
  if(arguments[arguments.length-1]=="#"){
     u=Url.split("#");
  }else{
     u=Url.split("?");
  }
  if (u.length==1){
    g='';
  }else{
    g=u[1];
  }
  if(g!=''){
     gg=g.split("&");
     var MaxI=gg.length;
     str = arguments[0]+"=";
     for(xm=0;xm<MaxI;xm++){
        if(gg[xm].indexOf(str)==0) {
          StrBack=gg[xm].replace(str,"");
          break;
        }
     }
  }
  return StrBack;
  }

function test(){
  this.x = 0;
  this.y = 0;
}
var uniFunction = {
    xMax : 9.0E+300,
    xMin : -9.0E+300,
    closedMax : false,
    closedMin : false,
    maxIter : 100,
    maxSearch : 25,
    xeps : 1.0E-006,
    feps : 1.0E-010,
    verbose : false
}


var BinomialAux = {
    xMax : 9.0E+300,
    xMin : -9.0E+300,
    closedMax : false,
    closedMin : false,
    maxIter : 100,
    maxSearch : 25,
    xeps : 1.0E-006,
    feps : 1.0E-010,
    verbose : false
}
BinomialAux.news = function(paramDouble1, paramDouble2){
    this.p0 = paramDouble1;
    this.n = paramDouble2;
    this.closedMin = (this.closedMax = true);
    this.xMin = 0.0;
    this.xMax = 1.0;
  }
BinomialAux.of = function(paramDouble){
   paramDouble = Math.max(0.0, Math.min(1.0, paramDouble));
      var d1 = Math.max(Math.abs(paramDouble - this.p0), 0.01);
      var d2 = Math.max(Math.sqrt(paramDouble * (1.0 - paramDouble) / this.n), d1 / 1000.0);
      var x = (paramDouble - this.p0) / d2;
      return x;
}

var PowerCalculatorAux = {
    xMax : 9.0E+300,
    xMin : -9.0E+300,
    closedMax : false,
    closedMin : false,
    maxIter : 100,
    maxSearch : 25,
    xeps : 1.0E-006,
    feps : 1.0E-010,
    verbose : false,
    xName:'',
    yName:'',
    o:null
}

PowerCalculatorAux.news = function(df, risk, obj){
  this.o = obj;
  this.xName = df;
  this.yName = risk;
}

PowerCalculatorAux.of = function(paramDouble){
  return this.o.eval(this.o, this.yName, this.xName, paramDouble);
}

var Chi2Aux = {
  xMax : 9.0E+300,
  xMin : -9.0E+300,
  closedMax : false,
  closedMin : false,
  maxIter : 100,
  maxSearch : 25,
  xeps : 1.0E-006,
  feps : 1.0E-010,
  verbose : false,
  df:0.0,
  lambda:0.0  
}


Chi2Aux.news = function(paramDouble1, paramDouble2){
  this.df = paramDouble1;
  this.lambda = paramDouble2;
  this.xMin = 0.0;
  this.closedMin = true;
  this.xeps = 1.E-005;
}

Chi2Aux.of = function(paramDouble){
  return Chi2.cdf3(paramDouble, this.df, this.lambda);
}

var Chi2Aux2 = {
    xMax : 9.0E+300,
    xMin : -9.0E+300,
    closedMax : false,
    closedMin : false,
    maxIter : 100,
    maxSearch : 25,
    xeps : 1.0E-006,
    feps : 1.0E-010,
    verbose : false,
    df:0.0,
    alpha:0.0 
}


Chi2Aux2.news = function(paramDouble1, paramDouble2){
  this.df = paramDouble1;
  this.alpha = paramDouble2;
  this.xMin = 0.0;
  this.closedMin = true;
  this.xeps = 1.E-005;
}

Chi2Aux2.of = function(paramDouble){
  return Chi2.power3(paramDouble, this.df, this.alpha);
}
//end of Unifunction.js

//begin of Solve.js
var Solve = {}
function NaNs(){return 0.0/0.0;}
Solve.illinois = function(paramUniFunction, paramDouble1, paramDouble2, paramDouble3){
    if ((paramDouble2 + paramDouble3) == null) {
      return NaNs();
    }
    var d1 = paramUniFunction.of(paramDouble2) - paramDouble1;
    var d2 = paramUniFunction.of(paramDouble3) - paramDouble1;
    var x = this.illinois(paramUniFunction, paramDouble1, paramDouble2, paramDouble3, d1, d2);
    return x;
  }
  
Solve.illinois = function(paramUniFunction, paramDouble1, paramDouble2, paramDouble3, paramDouble4, paramDouble5){
    var i = 0;
    if ((paramDouble4 + paramDouble5) == null) {
      return NaNs();
    }
    if (paramDouble4 * paramDouble5 > 0.0) {
      return NaNs();
    }while (i++ < paramUniFunction.maxIter){
      if (Math.abs(paramDouble4 - paramDouble5) <= 1.0E-099) {
        return NaNs();
      }
      var d1 = paramDouble3 - paramDouble5 * (paramDouble3 - paramDouble2) / (paramDouble5 - paramDouble4);
      var d2 = paramUniFunction.of(d1) - paramDouble1;
      if (d2==null) {
        return NaNs();
      }
      if (d2 * paramDouble5 > 0.0){
        paramDouble4 /= 2.0;
      }else{
        paramDouble2 = paramDouble3;
        paramDouble4 = paramDouble5;
      }
      paramDouble3 = d1;
      paramDouble5 = d2;
      if ((Math.abs(paramDouble2 - paramDouble3) <= paramUniFunction.xeps) || (Math.abs(d2) <= paramUniFunction.feps)) {
        return d1;
      }
    }
    return NaNs();
  }
  
Solve.search = function(paramUniFunction, paramDouble1, paramDouble2, paramDouble3){
    var d1 = paramDouble2 + paramDouble3;
    var i = 0;
    if ((paramDouble2 + paramDouble1 + paramDouble3)==null) {
      return NaNs();
    }
    var d2 = paramUniFunction.of(paramDouble2) - paramDouble1;
    var d3 = paramUniFunction.of(d1) - paramDouble1;
    if (d3==null) {
      return NaNs();
    }
    while (d2 * d3 > 0.0){
      i++;
      if (i > paramUniFunction.maxSearch) {
        return NaNs();
      }
      paramDouble3 = -1.5 * (d1 - paramDouble2) * d3 / (d3 - d2);
      
      paramDouble2 = d1;
      d2 = d3;
      d1 += paramDouble3;
      if (d1 <= paramUniFunction.xMin) {
        d1 = paramUniFunction.closedMin ? paramUniFunction.xMin : paramUniFunction.xMin + 0.1 * (paramDouble2 - paramUniFunction.xMin);
      }
      if (d1 >= paramUniFunction.xMax) {
        d1 = paramUniFunction.closedMax ? paramUniFunction.xMax : paramUniFunction.xMax - 0.1 * (paramUniFunction.xMax - paramDouble2);
      }
      d3 = paramUniFunction.of(d1) - paramDouble1;
      if (d3==null) {
        return NaNs();
      }
    }
    var x = this.illinois(paramUniFunction, paramDouble1, paramDouble2, d1, d2, d3);
    return x;
  }
//end of Solve.js

//begin of pc.js
var pc = {
    javaVersion : 1.0,
      panels : new Array(),
      listeners : new Array(),
      actionSource : "init",
      sourceIndex : -1,
      actors:new Array()
}

/*pc.saveVars= function(){
  var arrayOfDouble = new Array();
  for (var i = 0; i < this.actors.length; i++)
  {
    PiComponent localPiComponent = (PiComponent)this.actors.elementAt(i);
    arrayOfDouble[i] = ((localPiComponent instanceof DoubleComponent) ? getDVar(localPiComponent.getName()) : getIVar(localPiComponent.getName()));
  }
  return arrayOfDouble;
}*/

pc.solve = function(paramPowerCalculatorAux, paramDouble1, paramDouble2, paramDouble3){
   //var arrayOfDouble = saveVars();
    var d = Solve.search(paramPowerCalculatorAux, paramDouble1, paramDouble2, paramDouble3);
   //restoreVars(arrayOfDouble);
    return d;
  }

/*pc.getDVar = function(risk){
    if (risk.endsWith("]")) {
        return getDVar(parseArray(risk));
      }
    Field localField = getClass().getField(risk);
   double n= localField.getDouble(this);
   return n;
}*/
//end of pc.js

//begin of T.js
var T = {
  saveTail:0,
  saveDelta:0.0,
  saveDf:0.0,
  saveMu:0.0,
  saveTol:0.0,
  saveSE:0.0
  }
function NaNs(){return 0.0/0.0;};
T.cdf = function(paramDouble1, paramDouble2, paramDouble3){
  var d1 = 0.398942280401433;
  var d2 = 5.E-009;
  var i = 400;
  if ((paramDouble1 + paramDouble2 + paramDouble3)==null) {
    return (0.0 / 0.0);
  }
  if (paramDouble2 <= 0.0) {
    return NaNs();
  }
  if (paramDouble1 < 0.0) {
    return 1.0 - this.cdf(-paramDouble1, paramDouble2, -paramDouble3);
  }
  var d10 = paramDouble3 * paramDouble3;
  var d4 = Math.exp(-d10 / 2.0) / 2.0;
  var d5 = d1 * 2.0 * d4 * paramDouble3;
  var d6 = 0.5 - d4;
  var d7 = paramDouble1 * paramDouble1;
  d7 /= (d7 + paramDouble2);
  var d8 = 0.5;
  var d9 = d8 * paramDouble2;
  var d13 = Beta.cdf3(d7, d8, d9);
  var d11 = 2.0 * MoreMath.beta(d8, d9) * Math.pow(d7, d8) * Math.pow(1.0 - d7, d9);
  var d12;
  var d14;
  if (Math.abs(paramDouble3) > 1.0E-008)  {
    d12 = Math.pow(1.0 - d7, d9);
    d14 = 1.0 - d12;
    d12 = d9 * d7 * d12;
  }else{
    d14 = d12 = 0.0;
  }
  var j = 1;
  var d3 = d4 * d13 + d5 * d14;
  while ((d6 * (d13 - d11) > d2) && (j <= i)){
    d8 += 1.0;
    d13 -= d11;
    d14 -= d12;
    d11 *= d7 * (d8 + d9 - 1.0) / d8;
    d12 *= d7 * (d8 + d9 - 0.5) / (d8 + 0.5);
    d4 *= d10 / (2 * j);
    d5 *= d10 / (2 * j + 1);
    d6 -= d4;
    j++;
    d3 += d4 * d13 + d5 * d14;
  }
  if (j > i) {
    return NaNs();
  }
  return d3 + Normal.cdf(-paramDouble3);
}

T.cdf2 = function(paramDouble1, paramDouble2){
  return this.cdf(paramDouble1, paramDouble2, 0.0);
}
var TAux = {
    xMax : 9.0E+300,
    xMin : -9.0E+300,
    closedMax : false,
    closedMin : false,
    maxIter : 100,
    maxSearch : 25,
    xeps : 1.0E-006,
    feps : 1.0E-010,
    verbose : false,
    df : 0.0,
    delta : 0.0
};
TAux.news = function(paramDouble1, paramDouble2){
  this.df = paramDouble1;
  this.delta = paramDouble2;
}

TAux.of = function(paramDouble){
  return T.cdf(paramDouble, this.df, this.delta);
}
var TAux2 = {
    xMax : 9.0E+300,
    xMin : -9.0E+300,
    closedMax : false,
    closedMin : false,
    maxIter : 100,
    maxSearch : 25,
    xeps : 1.0E-006,
    feps : 1.0E-010,
    verbose : false,
    df : 0.0,
    critval : 0.0,
    tail : 0
};
TAux2.news = function(paramInt, paramDouble1, paramDouble2){
    this.tail = paramInt;
    this.df = paramDouble2;
    if (paramInt > 0) {
      this.critval = T.quantile2(1.0 - paramDouble1, paramDouble2);
    } else {
      this.critval = T.quantile2(1.0 - 0.5 * paramDouble1, paramDouble2);
    }
  }

TAux2.of = function(paramDouble){
  if (this.tail > 0) {
        return 1.0 - T.cdf(this.critval, this.df, paramDouble);
      }
      return T.cdf(-this.critval, this.df, paramDouble) + 1.0 - T.cdf(this.critval, this.df, paramDouble);
}


T.quantile = function(paramDouble1, paramDouble2, paramDouble3){
  if ((paramDouble1 + paramDouble2 + paramDouble3)==null) {
    return (0.0 / 0.0);
  }
  if ((paramDouble1 <= 0.0) || (paramDouble1 >= 1.0)) {
    return NaNs();
  }
  TAux.news(paramDouble2, paramDouble3);
  var d = paramDouble3 + 5.0 * (Math.pow(paramDouble1, 0.14) - Math.pow(1.0 - paramDouble1, 0.14));
  return Solve.search(TAux, paramDouble1, d - 0.005, 0.01);
}

T.quantile2 = function(paramDouble1, paramDouble2){
  return this.quantile(paramDouble1, paramDouble2, 0.0);
}

T.power = function(paramDouble1, paramDouble2, paramInt, paramDouble3){
  if ((paramDouble1 + paramDouble2 + paramDouble3)==null) {
    return (0.0 / 0.0);
  }
  if (paramInt < 0) {
    return this.power(-paramDouble1, paramDouble2, 1, paramDouble3);
  }
  if (paramInt == 0)  {
    var d = this.quantile2(paramDouble3 / 2.0, paramDouble2);
    return this.cdf(d, paramDouble2, paramDouble1) + 1.0 - this.cdf(-d, paramDouble2, paramDouble1);
  }
  var d = -this.quantile2(paramDouble3, paramDouble2);
  return 1.0 - this.cdf(d, paramDouble2, paramDouble1);
}

T.delta = function(paramDouble1, paramDouble2, paramInt, paramDouble3){
  if ((paramDouble1 + paramDouble2 + paramDouble3)==null) {
    return (0.0 / 0.0);
  }
  if ((paramDouble1 <= 0.0) || (paramDouble1 >= 1.0)) {
    return NaNs();
  }
  if ((paramDouble3 <= 0.0) || (paramDouble3 >= 1.0)) {
    return NaNs();
  }
  if (paramInt < 0) {
    return -this.delta(paramDouble1, paramDouble2, 1, paramDouble3);
  }
  TAux2.news(paramInt, paramDouble3, paramDouble2);
  var d;
  if (paramInt > 0) {
    d = this.quantile2(1.0 - paramDouble3, paramDouble2) + this.quantile2(1.0 - paramDouble1, paramDouble2);
  } else {
    d = this.quantile2(1.0 - 0.5 * paramDouble3, paramDouble2) + this.quantile2(1.0 - paramDouble1, paramDouble2);
  }
  return Solve.search(TAux2, paramDouble1, d - 0.05, 0.1);
}

T.rocArea = function(paramDouble1, paramDouble2, paramInt, paramDouble3){
  this.saveDelta = paramDouble1;
  this.saveDf = paramDouble2;
  this.saveTail = paramInt;
  return NumAnal.integral(T.class, "power", 0.0, 1.0, paramDouble3, false, 0.0, 1.0);
}

T.rocArea = function(paramDouble1, paramDouble2, paramInt){
  return this.rocArea(paramDouble1, paramDouble2, paramInt, 0.0001);
}

T.power1 = function(paramDouble){
  return this.power(this.saveDelta, this.saveDf, this.saveTail, paramDouble);
}

T.powerEquiv5 = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4, paramDouble5){
  if (paramDouble2 <= 0.0) {
    return 0.0;
  }
  var d1 = this.quantile2(1.0 - paramDouble5, paramDouble4);
  var d2 = paramDouble1 / paramDouble3;
  var d3 = paramDouble2 / paramDouble3;
  

  var d6 = 0.5 * paramDouble4;
  var d7 = 0.0, d4, d8, d9;
  if (Math.abs(d1) < 0.0001) {
    return Normal.cdf(d3 - d2) - Normal.cdf(-d3 - d2);
  }
  if (d1 > 0.0) {
    d4 = d3 / d1;
  } else {
    d4 = (1.0 / 0.0);
  }
  d8 = Math.sqrt(Chi2.quantile3(0.9999, paramDouble4, 0.0) / paramDouble4);
  d4 = Math.min(d4, d8);
  var d5 = d4 / 100.0;
  for ( d9 = d5; d9 < d4; d9 += d5){
    var d10 = Normal.cdf(d3 - d2 - d1 * d9) - Normal.cdf(d1 * d9 - d3 - d2);
    
    var d11 = Math.pow(d9, paramDouble4 - 1.0) * Math.exp(-d6 * d9 * d9);
    d7 += d10 * d11;
  }
  d9 = d6 * Math.log(d6) - MoreMath.logGamma(d6);
  return d7 * 2.0 * d5 * Math.exp(d9);
}

T.powerEquiv = function(paramDouble){
  return this.powerEquiv5(this.saveMu, this.saveTol, this.saveSE, this.saveDf, paramDouble);
}

T.rocEquiv = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4, paramDouble5){
  this.saveMu = paramDouble1;
  this.saveTol = paramDouble2;
  this.saveSE = paramDouble3;
  this.saveDf = paramDouble4;
  return NumAnal.integral(T.class, "powerEquiv", 0.0, 1.0, paramDouble5, false, 0.0, 1.0);
}

T.rocEquiv = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4){
  return this.rocEquiv(paramDouble1, paramDouble2, paramDouble3, paramDouble4, 0.0001);
}
//end of T.js

//begin of Beta.js
var Beta = {}

Beta.cdf3 = function(paramDouble1, paramDouble2, paramDouble3){
    var d11 = 1.0E-008;
    
    var j = 500;
    var k = 0;var m = 0;var n = 0;
    if ((paramDouble2 <= 0.0) || (paramDouble3 <= 0.0)){
      return (0.0 / 0.0);
    }
    paramDouble1 = paramDouble1 > 1.0 ? 1.0 : paramDouble1 < 0.0 ? 0.0 : paramDouble1;
    if ((paramDouble1 == 0.0) || (paramDouble1 == 1.0)) {
      return paramDouble1;
    }
    var d7 = Math.pow(paramDouble1, paramDouble2) * Math.pow(1.0 - paramDouble1, paramDouble3) * MoreMath.beta(paramDouble2, paramDouble3);
    if (paramDouble2 < 1.5){
      m = 1;
      d7 *= paramDouble1 * (paramDouble2 + paramDouble3) / paramDouble2;
      paramDouble2 += 1.0;
    }
    if (paramDouble3 < 1.5){
      n = 1;
      d7 *= (1.0 - paramDouble1) * (paramDouble2 + paramDouble3) / paramDouble3;
      paramDouble3 += 1.0;
    }
    var d1;
    if (paramDouble1 >= (paramDouble2 - 1.0) / (paramDouble2 + paramDouble3 - 2.0)){
      k = 1;
      paramDouble1 = 1.0 - paramDouble1;
      d1 = paramDouble2;paramDouble2 = paramDouble3;paramDouble3 = d1;
    }
    var d2 = 1.0 / (1.0 - paramDouble1 * (paramDouble2 + paramDouble3) / (paramDouble2 + 1.0));
    var d8 = d2;
    var d10 = d2;
    var i = 1;
    var d3;
    do{
      d3 = d2;
      var d6 = paramDouble2 + 2 * i;
      var d4 = i * (paramDouble3 - i) * paramDouble1 / (d6 * (d6 - 1.0));
      var d5 = -(paramDouble2 + i) * (paramDouble2 + paramDouble3 + i) * paramDouble1 / (d6 * (d6 + 1.0));
      d8 = d2 + d4 * d8;
      d2 = d8 + d5 * d2;
      d10 = 1.0 + d4 * d10;
      var d9 = d10 + d5;
      d8 /= d9;
      d10 /= d9;
      d2 /= d9;
      i++;
    } while ((Math.abs(d2 - d3) >= d11 * d2) && (i <= j));
    d2 = d7 * d2 / paramDouble2;
    if (k != 0){
      d1 = paramDouble2;paramDouble2 = paramDouble3;paramDouble3 = d1;
      paramDouble1 = 1.0 - paramDouble1;
      d2 = 1.0 - d2;
    }
    if (n != 0){
      paramDouble3 -= 1.0;
      d7 = paramDouble3 * d7 / ((1.0 - paramDouble1) * (paramDouble2 + paramDouble3));
      d2 -= d7 / paramDouble3;
    }
    if (m != 0){
      paramDouble2 -= 1.0;
      d2 += d7 / (paramDouble1 * (paramDouble2 + paramDouble3));
    }
    return d2;
  }
  
Beta.oldCdf = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4){
    var d6 = 1.0E-008;
    
    var j = 500;
    if ((paramDouble2 <= 0.0) || (paramDouble3 <= 0.0) || (paramDouble4 < 0.0)){
      return (0.0 / 0.0);
    }
    var d1 = this.cdf3(paramDouble1, paramDouble2, paramDouble3);
    var d5;
    if (paramDouble4 == 0.0){
      d5 = d1;
    }else {
      var d2 = Math.pow(paramDouble1, paramDouble2) * Math.pow(1.0 - paramDouble1, paramDouble3) * MoreMath.beta(paramDouble2, paramDouble3) / paramDouble2;
      paramDouble4 /= 2.0;
      var d3 = Math.exp(-paramDouble4);
      var d4 = 1.0 - d3;
      d5 = d3 * d1;
      var i = 0;
      do{
        i++;
        d1 -= d2;
        d2 *= paramDouble1 * (paramDouble2 + paramDouble3 + i - 1.0) / (paramDouble2 + i);
        d3 *= paramDouble4 / i;
        d4 -= d3;
        d5 += d3 * d1;
      } while ((d4 * (d1 - d2) >= d6) && (i <= j));
    }
    return d5;
  }
  
Beta.cdf4 = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4){
    var j = 0;
    if ((paramDouble2 <= 0.0) || (paramDouble3 <= 0.0) || (paramDouble4 < 0.0)){
      return (0.0 / 0.0);
    }
    paramDouble4 /= 2.0;
    if (paramDouble4 < 1.0E-008) {
      return this.cdf3(paramDouble1, paramDouble2, paramDouble3);
    }
    if (paramDouble4 > 15.0){
      j = 1 + Chi2.quantile2(1.0E-008, 2.0 * paramDouble4) / 2;
      while ((j > 0) && (Poisson.cdf(j - 1, paramDouble4) > 1.0E-008)) {
        j--;
      }
    }
    var k = Chi2.quantile2(0.99999999, 2.0 * paramDouble4) / 2;
    while (Poisson.cdf(k, paramDouble4) < 0.99999999) {
      k++;
    }
    var i = k;
    var d1 = this.cdf3(paramDouble1, paramDouble2 + i, paramDouble3);
    var d2 = MoreMath.logGamma(paramDouble2 + paramDouble3 + i - 1.0) - MoreMath.logGamma(paramDouble2 + i - 1.0) - MoreMath.logGamma(paramDouble3) + (paramDouble2 + i - 1.0) * Math.log(paramDouble1) + paramDouble3 * Math.log(1.0 - paramDouble1);
    
    d2 = Math.exp(d2) / (paramDouble2 + i - 1.0);
    var d3 = Math.exp(-paramDouble4 + i * Math.log(paramDouble4) - MoreMath.logGamma(i + 1));
    
    var d4 = d3 * d1;
    d3 *= i / paramDouble4;
    for (i = k - 1; i >= j; i--)
    {
      d1 += d2;
      d2 *= (paramDouble2 + i) / paramDouble1 / (paramDouble2 + i + paramDouble3 - 1.0);
      d4 += d3 * d1;
      d3 *= i / paramDouble4;
    }
    return d4;
  }
  
var BetaAux = {
  xMax:9.0E+300,
  xMin:-9.0E+300,
  closedMax:false,
  closedMin:false,
  maxIter:100,
  maxSearch:25,
  xeps:1.0E-006,
  feps:1.0E-010,
  verbose:false,
 a:0.0,
 b:0.0,
 lambda:0.0
}

BetaAux.news = function(paramDouble1, paramDouble2, paramDouble3){
  this.a = paramDouble1;
  this.b = paramDouble2;
  this.lambda = paramDouble3;
  this.xMin = 0.0;this.closedMin = true;
  this.xMax = 1.0;this.closedMax = true;
}

BetaAux.of = function(paramDouble){
  var x = Beta.cdf4(paramDouble, this.a, this.b, this.lambda);
  return x;
}

Beta.quantile = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4){
    if (paramDouble1 * (1.0 - paramDouble1) == 0.0) {
      return paramDouble1;
    }
    var d1 = 4.91 * (Math.pow(paramDouble1, 0.14) - Math.pow(1.0 - paramDouble1, 0.14));
    
    var d2 = paramDouble2 + paramDouble4 / 2.0;
    
    var d3 = d2 / (d2 + paramDouble3) + d1 * Math.sqrt(d2 * paramDouble3 / Math.pow(d2 + paramDouble3, 3.0));
    d3 = Math.min(0.99, Math.max(0.01, d3));
    
    BetaAux.news(paramDouble2, paramDouble3, paramDouble4);
    var x = Solve.search(BetaAux, paramDouble1, d3, 0.01);
    return x;
  }
  
Beta.quantile3 = function(paramDouble1, paramDouble2, paramDouble3){
  var x = this.quantile(paramDouble1, paramDouble2, paramDouble3, 0.0);
    return x;
  }
//end of Beta.js

//begin of MoreMath.js
var MoreMath = {};

MoreMath.logGamma = function(paramDouble){
    if (paramDouble <= 0.0){
      return (0.0 / 0.0);
    }
    if (paramDouble < 15.0) {
      return this.logGamma(paramDouble + 1.0) - Math.log(paramDouble);
    }
    var d = 1.0 / (paramDouble * paramDouble);
    var x = (paramDouble - 0.5) * Math.log(paramDouble) - paramDouble + 0.9189385332046727 + (1.0 - d * (1.0 - d * (1.0 - 0.75 * d) / 3.5) / 30.0) / 12.0 / paramDouble;
    return x;
  }
  
MoreMath.gamma = function(paramDouble){
  var x = Math.exp(this.logGamma(paramDouble));
    return x;
  }
  
MoreMath.beta = function(paramDouble1, paramDouble2){
  var x = Math.exp(this.logGamma(paramDouble1 + paramDouble2) - this.logGamma(paramDouble1) - this.logGamma(paramDouble2));
    return x;
  }
//end of MoreMath.js

//begin of Binomial.js
 var Binomial = {};

 Binomial.cdf = function(paramInt1, paramInt2, paramDouble){
    if ((paramDouble < 0.0) || (paramDouble > 1.0)){
      return 2147483647.0;
    }
    if (paramInt1 < 0) {
      return 0.0;
    }
    if (paramInt1 >= paramInt2) {
      return 1.0;
    }
    var x = Beta.cdf3(1.0 - paramDouble, 0.0 + paramInt2 - paramInt1, paramInt1 + 1.0);
    return x;
  }
  
 Binomial.quantile = function(paramDouble1, paramInt, paramDouble2){
    if ((paramDouble1 < 0.0) || (paramDouble1 > 1.0))
    {
      return 2147483647;
    }
    if ((paramDouble2 < 0.0) || (paramDouble2 > 1.0))
    {
      return 2147483647;
    }
    var i = (paramInt * paramDouble2 + Normal.quantile(paramDouble1) * Math.sqrt(paramInt * paramDouble2 * (1.0 - paramDouble2)));
    i = Math.max(-1, Math.min(paramInt, i));
    var d = this.cdf3(i, paramInt, paramDouble2);
    if (d < paramDouble1)
    {
      do
      {
        i++;d = this.cdf(i, paramInt, paramDouble2);
      } while (d < paramDouble1 - 1.0E-006);
      i--;
    }
    else if (d > paramDouble1){
      do
      {
        i--;d = this.cdf(i, paramInt, paramDouble2);
      } while (d > paramDouble1 + 1.0E-006);
    }
    return i;
  }
  
 Binomial.power = function(paramDouble1, paramDouble2, paramInt1, paramInt2, paramDouble3){
    var i;
    var d1;
    var d2;
    if (paramInt2 < 0){
      i = Binomial.quantile(paramDouble3, paramInt1, paramDouble1);
      d1 = this.cdf(i, paramInt1, paramDouble1);
      d2 = this.cdf(i, paramInt1, paramDouble2);
    }else if (paramInt2 == 0){
      i = this.quantile(paramDouble3 / 2.0, paramInt1, paramDouble1);
      d1 = this.cdf(i, paramInt1, paramDouble1);
      d2 = this.cdf(i, paramInt1, paramDouble2);
      i = 1 + Binomial.quantile(1.0 - paramDouble3 / 2.0, paramInt1, paramDouble1);
      d1 += 1.0 - this.cdf(i, paramInt1, paramDouble1);
      d2 += 1.0 - this.cdf(i, paramInt1, paramDouble2);
    }else{
      i = 1 + Binomial.quantile(1.0 - paramDouble3, paramInt1, paramDouble1);
      d1 = 1.0 - this.cdf(i, paramInt1, paramDouble1);
      d2 = 1.0 - this.cdf(i, paramInt1, paramDouble2);
    }
    return new Array(d2, d1);
  }
  
 Binomial.waldPower = function(paramDouble1, paramDouble2, paramInt1, paramInt2, paramDouble3){
    var d1 = 0.0;var d2 = 0.0;
    var d6 = paramInt2 == 0 ? paramDouble3 / 2.0 : paramDouble3;
    

    var d7 = Normal.quantile(1.0 - d6);
    var d4;
    var d5;
    var d3;
    var i;
    if (paramInt2 <= 0)
    {
      BinomialAux.news(paramDouble1, 4.0 + paramInt1);
      d4 = Math.max(paramDouble1 - d7 * Math.sqrt(paramDouble1 * (1.0 - paramDouble1) / paramInt1), 0.0);
      d5 = Math.min(0.01, 1.0 - d4);
      d3 = Solve.search(BinomialAux, -d7, d4, d5);
      i = Math.floor((paramInt1 + 4) * d3 - 2.0);
      d2 = this.cdf(i, paramInt1, paramDouble2);
      d1 = this.cdf(i, paramInt1, paramDouble1);
    }
    if (paramInt2 >= 0)
    {
       BinomialAux.news(1.0 - paramDouble1, 4.0 + paramInt1);
      d4 = Math.max(1.0 - paramDouble1 - d7 * Math.sqrt(paramDouble1 * (1.0 - paramDouble1) / paramInt1), 0.0);
      d5 = Math.min(0.01, 1.0 - d4);
      d3 = Solve.search(BinomialAux, -d7, d4, d5);
      i = Math.floor((paramInt1 + 4) * d3 - 2.0);
      d2 += this.cdf(i, paramInt1, 1.0 - paramDouble2);
      d1 += this.cdf(i, paramInt1, 1.0 - paramDouble1);
    }
    return new Array(d2, d1);
  }
//end if Binomial.js

//begin of PowerCalculator.js
function PowerCalculator(){}

PowerCalculator.prototype.javaVersion = 1.0;
PowerCalculator.prototype.panels = new Array();
PowerCalculator.prototype.actors = new Array();
PowerCalculator.prototype.listeners = new Array();
PowerCalculator.prototype.actionSource = "init";
PowerCalculator.prototype.sourceIndex = -1;

function endsWith(str, suffix){
  if(str.substring(str.length-suffix.length, str.length) == suffix){
    return true;
  }
  return false;
}

PowerCalculator.prototype.eval = function(paramString1, paramString2, paramDouble){
  this.setVar_sd(paramString2, paramDouble);
  this.callMethodFor(paramString2);
  return this.getDVar_s(paramString1);
}

function parseArray(paramString){
  var i = paramString.indexOf("[");
  var j = paramString.indexOf("]");
  var str = paramString.substring(0, i);
  var localInteger = parseInt(paramString.substring(i + 1, j));
  return [ str, localInteger ];
}

PowerCalculator.prototype.setVar_ad = function(paramArrayOfObject, paramDouble){
  var str = paramArrayOfObject[0];
  var i = paramArrayOfObject[1];
  this[str][i]=paramDouble;
}

PowerCalculator.prototype.setVar_sd = function(paramString, paramDouble){
  if (endsWith(paramString, "]")){
    this.setVar_ad(parseArray(paramString), paramDouble);
    return;
  }
  this[paramString] = paramDouble;
}

PowerCalculator.prototype.callMethodFor = function(paramString){
  var localObject;
  if (endsWith(paramString, "]")){
    localObject = parseArray(paramString);
    this.actionSource = localObject[0];
    this.sourceIndex = localObject[1];
  }else{
    this.actionSource = paramString;
    this.sourceIndex = -1;
  }
  var x = this[this.actionSource + "_changed"];
  if(x!=undefined)
    x();
  this.click();
}

PowerCalculator.prototype.getDVar_a = function( paramArrayOfObject){
  var str = paramArrayOfObject[0];
  var i = paramArrayOfObject[1];
  var localObject = obj[str];
  if(localObject[i] != NaN)
    return localObject[i];
  return (0.0 / 0.0);
}

PowerCalculator.prototype.getDVar_s = function(risk){
  if (endsWith(risk, "]")) {
      return this.getDVar_a(parseArray(risk));
  }
  var n= obj[risk];
  if(n != NaN)
    return n;
  return (0.0 / 0.0);
}

PowerCalculator.prototype.getDVars_pddd = function(paramPowerCalculatorAux, paramDouble1, paramDouble2, paramDouble3){
  var arrayOfDouble = this.saveVars();
  var d = Solve.search(paramPowerCalculatorAux, paramDouble1, paramDouble2, paramDouble3);
  this.restoreVars(arrayOfDouble);
  return d;
}
PowerCalculator.prototype.saveVars = function(){
  var arrayOfDouble = new Array(this.actors.size());
//  for (var i = 0; i < this.actors.size(); i++){
//    var localPiComponent = this.actors[i];
//    arrayOfDouble[i] = ((localPiComponent instanceof DoubleComponent) ? getDVar(localPiComponent.getName()) : getIVar(localPiComponent.getName()));
//  }
  return arrayOfDouble;
}

PowerCalculator.prototype.restoreVars = function(paramArrayOfDouble){
//  for (var i = 0; i < this.actors.size(); i++)  {
//    PiComponent localPiComponent = (PiComponent)this.actors.elementAt(i);
//    if ((localPiComponent instanceof DoubleComponent)){
//      setVar(localPiComponent.getName(), paramArrayOfDouble[i]);
//      ((DoubleComponent)localPiComponent).setValue(paramArrayOfDouble[i]);
//    }
//    else
//    {
//      setVar(localPiComponent.getName(), (int)paramArrayOfDouble[i]);
//      ((IntComponent)localPiComponent).setValue((int)paramArrayOfDouble[i]);
//    }
//  }
}

PowerCalculator.prototype.solve = function(paramPowerCalculatorAux, paramDouble1, paramDouble2, paramDouble3){
  var arrayOfDouble = this.saveVars();
  var d = Solve.search(paramPowerCalculatorAux, paramDouble1, paramDouble2, paramDouble3);
  restoreVars(arrayOfDouble);
  return d;
}
//end of PowerCalculator.js

//begin of OnePGUI.js
function OnePGUI(){};

OnePGUI.prototype = new PowerCalculator();
OnePGUI.prototype.p0=0.0;
OnePGUI.prototype.p=0.0;
OnePGUI.prototype.n=0.0;
OnePGUI.prototype.alt=0;
OnePGUI.prototype.Alpha=0.0;
OnePGUI.prototype.Method=0;
OnePGUI.prototype.sizes=0.0;
OnePGUI.prototype.Power=0.0;

OnePGUI.prototype.gui = function(){
  this.p0=0.5;
  this.p=0.6;
  this.n=50.0;
  this.alt=1;
  this.Alpha=0.05;
  this.Method=1;
  this.sizes=0.06;
  this.Power=0.0;
}

OnePGUI.prototype.nPower = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4, paramInt, paramDouble5){
    var d1;
    if (paramInt > 0)
    {
      d1 = Normal.quantile3(1.0 - paramDouble5, paramDouble1, paramDouble2);
      return 1.0 - Normal.cdf3(d1, paramDouble3, paramDouble4);
    }
    if (paramInt < 0)
    {
      d1 = Normal.quantile3(paramDouble5, paramDouble1, paramDouble2);
      return Normal.cdf3(d1, paramDouble3, paramDouble4);
    }
    var d2 = Normal.quantile3(paramDouble5 / 2.0, paramDouble1, paramDouble2);
    var d3 = Normal.quantile3(1.0 - paramDouble5 / 2.0, paramDouble1, paramDouble2);
    var x = Normal.cdf3(d2, paramDouble3, paramDouble4) - Normal.cdf3(d3, paramDouble3, paramDouble4);
    x = 1.0+x;
    return x;
  }
OnePGUI.prototype.click = function(){
  var arrayOfDouble = new Array();
  var i = this.alt - 1;
      this.n = Math.max(2.0,Math.round(this.n));
      this.p0 = Math.min(Math.max(this.p0, 0.01), 0.99);
      this.p = Math.min(Math.max(this.p, 0.01), 0.99);
      switch (this.Method){
      case 0: 
        arrayOfDouble = Binomial.power(this.p0, this.p, this.n, i, this.Alpha);
        this.Power = arrayOfDouble[0];
        this.sizes = arrayOfDouble[1];
        break;
      case 1: 
        var d1 = Math.sqrt(this.p0 * (1.0 - this.p0) / this.n);
        var d2 = Math.sqrt(this.p * (1.0 - this.p) / this.n);
        var xx = this.nPower(this.p0, d1, this.p, d2, i, this.Alpha);
        this.Power = xx;
        this.sizes = this.Alpha;
        break;
      case 2: 
        this.Power = this.bPower((this.n - 1.0) * this.p0, (this.n - 1.0) * (1.0 - this.p0), (this.n - 1.0) * this.p, (this.n - 1.0) * (1.0 - this.p), i, this.Alpha);
        
        this.sizes = this.Alpha;
        break;
      case 3: 
        arrayOfDouble = Binomial.waldPower(this.p0, this.p, this.n, i, this.Alpha);
        this.Power = arrayOfDouble[0];
        this.sizes = arrayOfDouble[1];
      }
}

OnePGUI.prototype.worstCase_changed = function() {
  this.click();
}

OnePGUI.prototype.bPower = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4, paramInt, paramDouble5){
  var d1;
  if (paramInt > 0){
    d1 = Beta.quantile3(1.0 - paramDouble5, paramDouble1, paramDouble2);
    return 1.0 - Beta.cdf3(d1, paramDouble3, paramDouble4);
  }
  if (paramInt < 0)
  {
    d1 = Beta.quantile3(paramDouble5, paramDouble1, paramDouble2);
    return Beta.cdf3(d1, paramDouble3, paramDouble4);
  }
  var d2 = Beta.quantile3(paramDouble5 / 2.0, paramDouble1, paramDouble2);
  var d3 = Beta.quantile3(1.0 - paramDouble5 / 2.0, paramDouble1, paramDouble2);
  return 1.0 + Beta.cdf3(d2, paramDouble3, paramDouble4) - Beta.cdf3(d3, paramDouble3, paramDouble4);
}
OnePGUI.prototype.Power_changed = function(){
  this.Power = Math.min(0.99, Math.max(this.Alpha, this.Power));
  PowerCalculatorAux.news("n", "Power", this);
  PowerCalculatorAux.closedMin = true;
  PowerCalculatorAux.xMin = 2.0;
  PowerCalculatorAux.xeps = 0.5;
  
  this.n = pc.solve(PowerCalculatorAux, this.Power, this.n, 20.0);
  console.log(onePGUI.Power, onePGUI.n)
  this.click();
}

function StringNaN(o){
  if(o!=null && o!="" && o!=undefined)
    return false;
  return true;
}
var onePGUI = null;
module.exports.OnePGUI_cfap = function(id, key, value){
  if(isNullOrEmpty(onePGUI)){
    onePGUI = new OnePGUI();
    onePGUI.gui();
    onePGUI.click();
  }
  var mm = 0;
  if(id!=null && id!=''){
    if(key!=null && key!=""){
      if(key=="p0"){
          if(!StringNaN(value)){
          onePGUI.p0 = parseFloat(value);
          }
          onePGUI.click();
        }
        if(key=="p"){
          if(!StringNaN(value)){
          onePGUI.p = parseFloat(value);
          }
          onePGUI.click();
        }
        if(key=="ssize"){
          if(!StringNaN(value)){
          onePGUI.n = parseFloat(value);
          }
        onePGUI.click();
        }
        if(key=="altt"){
          if(!StringNaN(value)){
            onePGUI.alt= parseInt(value);
          }
        onePGUI.click();
        }
        if(key=="alpha"){
          if(!StringNaN(value)){
            onePGUI.Alpha = parseFloat(value);
          }
        onePGUI.click();
        }
        if(key=="method"){
          if(!StringNaN(value)){
            onePGUI.Method = parseInt(value);
          }
        onePGUI.click();
        }
      if(key=="power"){
        if(!StringNaN(value)){
          onePGUI.Power = parseFloat(value);
        }
        onePGUI.Power_changed();
      }
      }
}else{
  onePGUI = new OnePGUI();
  onePGUI.gui();
  onePGUI.click();
}
  return onePGUI;
}
//end of OnePGUI.js

//begin of OneTGUI.js
var OneTGUI = {
  sigma:0.0,
  n:0.0,
  diff:0.0,
  alpha:0.0,
  df:0.0,
  power:0.0,
  delta:0.0,
  eqs:0,
  eqn:0,
  tt:0,
  opt:0
}
OneTGUI.news = function(){
  this.sigma=0.0;
  this.n=0.0;
  this.diff=0.0;
  this.alpha=0.0;
  this.df=0.0;
  this.power=0.0;
  this.delta=0.0;
  this.eqs=0;
  this.eqn=0;
  this.tt=0;
  this.opt=0;
}
OneTGUI.gui = function(){
  this.sigma=1.0;
  this.diff=0.5;
  this.n=25.0;
  this.power=0.5; 
  this.opt=0;
  this.alpha=0.05;
  this.tt=1;
}

OneTGUI.click = function(){
  this.n = Math.max(Math.round(this.n), 2.0);
  this.delta = (Math.sqrt(this.n) * this.diff / this.sigma);
  this.power = T.power(this.delta, this.n - 1.0, 1 - this.tt, this.alpha);
}

OneTGUI.power_changed = function(){
  switch (this.opt){
  case 0: 
    if (Math.abs(this.diff) < 0.001 * this.sigma) {
      return;
    }
    for (var i = 0; i < 3; i++){
      var d2 = this.n;
      this.delta = T.delta(this.power, this.n - 1.0, 1 - this.tt, this.alpha);
      this.n = Math.pow(this.delta * this.sigma / this.diff, 2.0);
      if (isNaN(this.n)){
        this.n = d2;
        break;
      }
    }
    break;
  case 1: 
    var d1 = this.diff;
    this.delta = T.delta(this.power, this.n - 1.0, 1 - this.tt, this.alpha);
    this.diff = (this.sigma * this.delta / Math.sqrt(this.n));
    if (isNaN(this.diff)) {
      this.diff = d1;
    }
    break;
  }
  this.click();
}
module.exports.handle = function(id, key, value){
     if(isNullOrEmpty(OneTGUI.n)){
       OneTGUI.news();
       OneTGUI.gui();
       OneTGUI.click();
     }
     try{
     if(!isNullOrEmpty(id)){
       if(!isNullOrEmpty(key)){
         if(key=="sigma"){
           if(!isNullOrEmpty(value)){
             OneTGUI.sigma=parseFloat(value);
             OneTGUI.click();
           }
         }
         if(key=="diff"){
           if(!isNullOrEmpty(value)){
           OneTGUI.diff=parseFloat(value);
             OneTGUI.click();
           }
         }
         if(key=="n"){
           if(!isNullOrEmpty(value)){
             OneTGUI.n=parseFloat(value);
             OneTGUI.click();
           }
         }
         if(key=="power"){
           if(!isNullOrEmpty(value)){
             OneTGUI.power=parseFloat(value);
             OneTGUI.power_changed();
           }
         }
         if(key=="opt"){
           if(!isNullOrEmpty(value)){
             OneTGUI.opt=parseInt(value);
             OneTGUI.click();
           }
         }
         if(key=="alpha"){
           if(!isNullOrEmpty(value)){
             OneTGUI.alpha=parseFloat(value);
             OneTGUI.click();
           }
         }
         if(key=="tt"){
            if(!isNullOrEmpty(value)){
              OneTGUI.tt=parseInt(value);
            }else{
              OneTGUI.setTt(0);
            }
            OneTGUI.click();
         }
       }
     }
     }catch(e){
       OneTGUI.gui();
       OneTGUI.click();
     }
     return OneTGUI;
}
//end of OneTGUI.js

//begin of Chi2.js
var Chi2 = {
  saveDf : 1.0,
  saveLambda : 0.0
}

Chi2.power3 = function(paramDouble1, paramDouble2, paramDouble3){
  if ((paramDouble3 <= 0.0) || (paramDouble3 >= 1.0)) {
    return (0.0 / 0.0);
  }
  if (paramDouble1 < 0.0) {
    return (0.0 / 0.0);
  }
  if (paramDouble2 < 0.01)  {
    return (0.0 / 0.0);
  }
  var d = this.quantile3(1.0 - paramDouble3, paramDouble2, 0.0);
  return 1.0 - this.cdf3(d, paramDouble2, paramDouble1);
}

Chi2.quantile3 = function(paramDouble1, paramDouble2, paramDouble3){
  if ((paramDouble1 < 0.0) || (paramDouble1 >= 1.0))  {
    return (0.0 / 0.0);
  }
  if (paramDouble3 < 0.0)  {
    return (0.0 / 0.0);
  }
  if (paramDouble2 < 0.01)  {
    return (0.0 / 0.0);
  }
  var d = paramDouble2 + paramDouble3 + Normal.quantile(paramDouble1) * Math.sqrt(2.0 * paramDouble2 + 4.0 * paramDouble3);
  Chi2Aux.news(paramDouble2, paramDouble3);
  return Solve.search(Chi2Aux, paramDouble1, d, 0.1);
}
Chi2.quantile2 = function(paramDouble1, paramDouble2){
  return this.quantile3(paramDouble1, paramDouble2, 0.0);
}
Chi2.cdf3 = function(paramDouble1, paramDouble2, paramDouble3){
  var d12 = 1.0E-008;
  var j = 500;
  if (paramDouble1 <= 0.0) {
    return 0.0;
  }
  if (paramDouble2 < 0.01)  {
    return (0.0 / 0.0);
  }
  if (paramDouble3 < 0.0)  {
    return (0.0 / 0.0);
  }
  var d4 = 0.0;
  var d1 = paramDouble1 / 2.0;
  var d2 = paramDouble2 / 2.0;
  var d6 = Math.exp(d2 * Math.log(d1) - d1 - MoreMath.logGamma(d2));
  var i = 1;
  var d7;
  var d10;
  var d3;
  if (paramDouble1 / paramDouble2 > 1.0)  {
    d7 = 0.0;
    var d8 = d6 / d1;
    d10 = 1.0 / d1;
    do    {
      i++;
      var d9;
      var d11;
      if (2 * (i / 2) < i)
      {
        d9 = (i - 1.0) / 2.0;
        d11 = d1;
      }      else      {
        d9 = (i - paramDouble2) / 2.0;
        d11 = 1.0;
      }
      d10 = 1.0 / (d9 * d10 + d11);
      d3 = d8;
      d8 = d10 * (d9 * d7 + d11 * d8);
      d7 = d10 * d3;
      d3 = d4;
      d4 = 1.0 - d8;
    } while ((Math.abs(d3 - d4) >= d12) && (i <= j));
    if (i > j) {}
  } else {
    d3 = d6 / d2;
    d4 = d3;
    do  {
      d3 *= paramDouble1 / (paramDouble2 + 2 * i);
      d4 += d3;
      i++;
    } while ((d3 >= d12) && (i <= j));
    if (i > j) {  }
  }
  if (paramDouble3 > 0.0)  {
    var d5 = Math.exp(-paramDouble3 / 2.0);
    d10 = 1.0 - d5;
    d7 = d4 - d6 / d2;
    d6 = d6 * d1 / (d2 * (d2 + 1.0));
    d2 += 2.0;
    d4 *= d5;
    i = 0;
    do    {
      i++;
      d5 *= paramDouble3 / (2 * i);
      d10 -= d5;
      d4 += d5 * d7;
      d7 -= d6;
      d6 *= d1 / d2;
      d2 += 1.0;
    } while ((d10 * d7 >= d12) && (i <= j));
    if (i > j) {
    }
  }
  return d4;
}

Chi2.cdf2 = function(paramDouble1, paramDouble2){
  return this.cdf3(paramDouble1, paramDouble2, 0.0);
}

Chi2.lambda = function (paramDouble1, paramDouble2, paramDouble3){
  if ((paramDouble3 <= 0.0) || (paramDouble3 >= 1.0))  {
    return (0.0 / 0.0);
  }
  if ((paramDouble1 <= 0.0) || (paramDouble1 >= 1.0))  {
    return (0.0 / 0.0);
  }
  if (paramDouble2 < 0.01)  {
    return (0.0 / 0.0);
  }
  var d1 = Math.pow(Normal.quantile(1.0 - paramDouble1), 2.0);
  var d2 = this.quantile3(1.0 - paramDouble3, paramDouble2, 0.0);
  var d3 = 2.0 * (d2 - paramDouble2) + 4.0 * d1;
  var d4 = Math.pow(d2 - paramDouble2, 2.0) - 2.0 * paramDouble2 * d1;
  var d5 = Math.pow(d3, 2.0) - 4.0 * d4;
  var d6 = paramDouble1 > 0.5 ? 1.0 : -1.0;
  var d7 = d5 < 0.0 ? d2 : (d3 + d6 * Math.sqrt(d5)) / 2.0;
  d7 = Math.max(0.5, d7);
  Chi2Aux2.news(paramDouble2, paramDouble3);
  return Solve.search(Chi2Aux2, paramDouble1, d7, 0.5);
}

Chi2.rocArea3 = function (paramDouble1, paramDouble2, paramDouble3){
  saveLambda = paramDouble1;
  saveDf = paramDouble2;
  return NumAnal.integral(Chi2.class, "power", 0.0, 1.0, paramDouble3, false, 0.0, 1.0);
}

Chi2.rocArea2 = function (paramDouble1, paramDouble2)
{
  return this.rocArea3(paramDouble1, paramDouble2, 0.0001);
}

Chi2.power1 = function (paramDouble){
  return this.power3(this.saveLambda, this.saveDf, paramDouble);
}
//end of Chi2.js

//begin of Tukey.js
var Tukey = {
  M_1_SQRT_2PI : 0.398942280401433,
    M_LN2 : 0.6931471805599453
}
Tukey.wprob = function(paramDouble1, paramDouble2, paramDouble3){
    var i = 12;
    var j = 6;
    var d1 = 1.0;
    var d2 = -30.0;
    var d3 = -50.0;
    var d4 = 60.0;
    var d5 = 8.0;
    var d6 = 3.0;
    var d7 = 2.0;
    var d8 = 3.0;
    var arrayOfDouble1 = [ 0.9815606342467192, 0.9041172563704749, 0.7699026741943047, 0.5873179542866175, 0.3678314989981802, 0.1252334085114689 ];
    var arrayOfDouble2 = [ 0.04717533638651183, 0.1069393259953184, 0.1600783285433462, 0.2031674267230659, 0.2334925365383548, 0.2491470458134028 ];
    var d23 = paramDouble1 * 0.5;
    if (d23 >= d5) {
      return 1.0;
    }
    var d11 = 2.0 * Normal.cdf(d23) - 1.0;
    if (d11 >= Math.exp(d3 / paramDouble3)) {
      d11 = Math.pow(d11, paramDouble3);
    } else {
      d11 = 0.0;
    }
    var d26;
    if (paramDouble1 > d6) {
      d26 = d7;
    } else {
      d26 = d8;
    }
    var d14 = d23;
    var d13 = (d5 - d23) / d26;
    var d15 = d14 + d13;
    var d18 = 0.0;

    var d17 = paramDouble3 - 1.0;
    for (var d25 = 1.0; d25 <= d26; d25 += 1.0)    {
      var d19 = 0.0;
      var d9 = 0.5 * (d15 + d14);

      var d12 = 0.5 * (d15 - d14);
      for (var m = 1; m <= i; m++) {
        var k;
        var d27;
        if (j < m){
          k = i - m + 1;
          d27 = arrayOfDouble1[(k - 1)];
        } else{
          k = m;
          d27 = -arrayOfDouble1[(k - 1)];
        }
        var d16 = d12 * d27;
        var d10 = d9 + d16;

        var d22 = d10 * d10;
        if (d22 > d4) {
          break;
        }
        var d21 = 2.0 * Normal.cdf(d10);
        var d20 = 2.0 * Normal.cdf3(d10, paramDouble1, 1.0);

        var d24 = d21 * 0.5 - d20 * 0.5;
        if (d24 >= Math.exp(d2 / d17)) {
          d24 = arrayOfDouble2[(k - 1)] * Math.exp(-(0.5 * d22)) * Math.pow(d24, d17);
          d19 += d24;
        }
      }
      d19 *= 2.0 * d12 * paramDouble3 * 0.398942280401433;
      d18 += d19;
      d14 = d15;
      d15 += d13;
    }
    d11 = d18 + d11;
    if (d11 <= Math.exp(d2 / paramDouble2)) {
      return 0.0;
    }
    d11 = Math.pow(d11, paramDouble2);
    if (d11 >= d1) {
      d11 = 1.0;
    }
    return d11;
  }
  
Tukey.cdf4 = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4){
  var i = 16;
  var j = 8;
    var d1 = -30.0;
    var d2 = 1.0E-014;
    var d3 = 100.0;
    var d4 = 800.0;
    var d5 = 5000.0;
    var d6 = 25000.0;
    var d7 = 1.0;
    var d8 = 0.5;
    var d9 = 0.25;
    var d10 = 0.125;
    var arrayOfDouble1 = [ 0.9894009349916499, 0.944575023073233, 0.8656312023878318, 0.755404408355003, 0.6178762444026438, 0.4580167776572274, 0.2816035507792589, 0.09501250983763744 ];
    var arrayOfDouble2 = [ 0.0271524594117541, 0.06225352393864789, 0.09515851168249279, 0.1246289712555339, 0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.189450610455069 ];
    if ((isNullOrEmpty(paramDouble1)) || (isNullOrEmpty(paramDouble4)) || (isNullOrEmpty(paramDouble2)) || (isNullOrEmpty(paramDouble3))) {
      return (0.0 / 0.0);
    }
    if (paramDouble1 <= 0.0) {
      return 0.0;
    }
    if ((paramDouble3 < 2.0) || (paramDouble4 < 1.0) || (paramDouble2 < 2.0)) {
      return (0.0 / 0.0);
    }
    if (paramDouble1==Number.POSITIVE_INFINITY) {
      return 1.0;
    }
    if (paramDouble3 > d6) {
      return this.wprob(paramDouble1, paramDouble4, paramDouble2);
    }
    var d12 = paramDouble3 * 0.5;
    var d14 = d12 * Math.log(paramDouble3) - paramDouble3 * 0.6931471805599453 - MoreMath.logGamma(d12);
    var d13 = d12 - 1.0;
    var d15 = paramDouble3 * 0.25;
    var d21;
    if (paramDouble3 <= d3) {
      d21 = d7;
    } else if (paramDouble3 <= d4) {
      d21 = d8;
    } else if (paramDouble3 <= d5) {
      d21 = d9;
    } else {
      d21 = d10;
    }
    d14 += Math.log(d21);
    var d16;
    var d11 = d16 = 0.0;
    for (var k = 1; k <= 50; k++)    {
      d16 = 0.0;
      var d20 = (2 * k - 1) * d21;
      for (var n = 1; n <= i; n++){
        var m;
        var d19;
        if (j < n){
          m = n - j - 1;
          d19 = d14 + d13 * Math.log(d20 + arrayOfDouble1[m] * d21) - (arrayOfDouble1[m] * d21 + d20) * d15;
        }else{
          m = n - 1;
          d19 = d14 + d13 * Math.log(d20 - arrayOfDouble1[m] * d21) + (arrayOfDouble1[m] * d21 - d20) * d15;
        }
        if (d19 >= d1){
          var d17;
          if (j < n) {
            d17 = paramDouble1 * Math.sqrt((arrayOfDouble1[m] * d21 + d20) * 0.5);
          } else {
            d17 = paramDouble1 * Math.sqrt((-(arrayOfDouble1[m] * d21) + d20) * 0.5);
          }
          var d22 = wprob(d17, paramDouble4, paramDouble2);
          var d18 = d22 * arrayOfDouble2[m] * Math.exp(d19);
          d16 += d18;
        }
      }
      if ((k * d21 >= 1.0) && (d16 <= d2)) {
        break;
      }
      d11 += d16;
    }
    if (d16 > d2) {
      return (0.0/0.0);
    }
    if (d11 > 1.0) {
      d11 = 1.0;
    }
    return d11;
  }
  
Tukey.cdf3 = function(paramDouble1, paramDouble2, paramDouble3) {
    return this.cdf4(paramDouble1, paramDouble2, paramDouble3, 1.0);
  }
  
Tukey.qinv = function(paramDouble1, paramDouble2, paramDouble3){
    var d1 = 0.322232421088;
    var d2 = 0.099348462606;
    var d3 = -1.0;
    var d4 = 0.588581570495;
    var d5 = -0.342242088547;
    var d6 = 0.531103462366;
    var d7 = -0.204231210125;
    var d8 = 0.10353775285;
    var d9 = -4.53642210148E-005;
    var d10 = 0.0038560700634;
    var d11 = 0.8832;
    var d12 = 0.2368;
    var d13 = 1.214;
    var d14 = 1.208;
    var d15 = 1.4142;
    var d16 = 120.0;
    var d17 = 0.5 - 0.5 * paramDouble1;
    var d20 = Math.sqrt(Math.log(1.0 / (d17 * d17)));
    var d19 = d20 + ((((d20 * d9 + d7) * d20 + d5) * d20 + d3) * d20 + d1) / ((((d20 * d10 + d8) * d20 + d6) * d20 + d4) * d20 + d2);
    if (paramDouble3 < d16) {
      d19 += (d19 * d19 * d19 + d19) / paramDouble3 / 4.0;
    }
    var d18 = d11 - d12 * d19;
    if (paramDouble3 < d16) {
      d18 += -d13 / paramDouble3 + d14 * d19 / paramDouble3;
    }
    return d19 * (d18 * Math.log(paramDouble2 - 1.0) + d15);
  }
  
Tukey.quantile4 = function(paramDouble1, paramDouble2, paramDouble3, paramDouble4){
  var d1 = 0.0001;
    var i = 50;
    var d2 = 0.0;
    if ((isNullOrEmpty(paramDouble1)) || (isNullOrEmpty(paramDouble4)) || (isNullOrEmpty(paramDouble2)) || (isNullOrEmpty(paramDouble3))) {
      return (0.0 / 0.0);
    }
    if ((paramDouble1 < 0.0) || (paramDouble1 >= 1.0)) {
      return (0.0 / 0.0);
    }
    if ((paramDouble3 < 2.0) || (paramDouble4 < 1.0) || (paramDouble2 < 2.0)) {
      return (0.0 / 0.0);
    }
    if (paramDouble1 == 0.0) {
      return 0.0;
    }
    var d5 = this.qinv(paramDouble1, paramDouble2, paramDouble3);
    var d3 = this.cdf4(d5, paramDouble2, paramDouble3, paramDouble4) - paramDouble1;
    var d6;
    if (d3 > 0.0) {
      d6 = Math.max(0.0, d5 - 1.0);
    } else {
      d6 = d5 + 1.0;
    }
    var d4 = this.cdf4(d6, paramDouble2, paramDouble3, paramDouble4) - paramDouble1;
    for (var j = 1; j < i; j++)
    {
      d2 = d6 - d4 * (d6 - d5) / (d4 - d3);
      d3 = d4;
      d5 = d6;
      if (d2 < 0.0){
        d2 = 0.0;
        d4 = -paramDouble1;
      }
      d4 = this.cdf4(d2, paramDouble2, paramDouble3, paramDouble4) - paramDouble1;
      d6 = d2;
      var d7 = Math.abs(d6 - d5);
      if (d7 < d1) {
        return d2;
      }
    }
    return d2;
  }
  
Tukey.quantile3 = function (paramDouble1, paramDouble2, paramDouble3)  {
    return this.quantile4(paramDouble1, paramDouble2, paramDouble3, 1.0);
  }
//end of Tukey.js

//begin of Poisson.js
var Poisson = {}

Poisson.cdf = function(paramInt, paramDouble)  {
    if (paramInt < 0) {
      return 0.0;
    }
    if (paramDouble <= 0.0)
    {
      return (0.0 / 0.0);
    }
    return 1.0 - Chi2.cdf2(2.0 * paramDouble, 2.0 * (paramInt + 1.0));
  }
  
Poisson.quantile = function (paramDouble1, paramDouble2)  {
    if (paramDouble2 <= 0.0)    {
      return 2147483647;
    }
    if ((paramDouble1 < 0.0) || (paramDouble1 >= 1.0))    {
      return 2147483647;
    }
    if (paramDouble1 == 0.0) {
      return -1;
    }
    var i = (paramDouble2 + Normal.quantile(paramDouble1) * Math.sqrt(paramDouble2));
    i = Math.max(-1, i);
    var d = this.cdf(i, paramDouble2);
    if (d < paramDouble1)    {
      do      {
        i++;d = this.cdf(i, paramDouble2);
      } while (d < paramDouble1 - 1.0E-006);
      i--;
    }else if (d > paramDouble1){
      do {
        i--;d = this.cdf(i, paramDouble2);
      } while (d > paramDouble1 + 1.0E-006);
    }
    return i;
  }

//end of Poisson.js

//begin of Pilot.js
function Pilot(){}
  
Pilot.prototype = new PowerCalculator();
Pilot.prototype.risk = 0.0;
Pilot.prototype.pctUnder = 0.0;
Pilot.prototype.df = 0.0;

Pilot.prototype.gui = function (){
      this.pctUnder=20.0;
      this.risk=0.1;
      this.df=80.0;
    }
    
Pilot.prototype.afterShow = function() {    }
    
   
Pilot.prototype.click = function (){
      this.df =Math.round(this.df);
      this.risk = Chi2.cdf2((1.0 - 0.01 * this.pctUnder) * this.df, this.df);
    }
    
Pilot.prototype.risk_changed = function (){
    this.df = pilot.solve( "df","risk", this.risk, this.df, 0.1 * Math.max(10.0, this.df) );
      this.click();
    }
//
Pilot.prototype.solve = function(df, risk, paramDouble1, paramDouble2, paramDouble3){
  PowerCalculatorAux.news(df, risk, this);
  var d = Solve.search(PowerCalculatorAux, paramDouble1, paramDouble2, paramDouble3);
  return d;
}

var pilot = null;
module.exports.pilot_handle = function(id, key, value){
  if(isNullOrEmpty(pilot)){
    pilot = new Pilot();
    pilot.gui();
    pilot.click();
    }
    if(!isNullOrEmpty(id)){
          if(!isNullOrEmpty(key)){
            if(key=="pctUnder"){
              if(!isNullOrEmpty(value)){
                pilot.pctUnder = parseFloat(value);
                pilot.click();
              }
            }
            if(key=="risk"){
              if(!isNullOrEmpty(value)){
                pilot.risk = parseFloat(value);
                pilot.risk_changed();
               
              }
            }
            if(key=="df"){
              if(!isNullOrEmpty(value)){
                pilot.ff = parseFloat(value);
                pilot.click();
              }
            }
          }
    }else{
      pilot = new Pilot();
      pilot.gui();
    }
    return pilot;
}
//End if Pilot.js

//begin of Rsquare.js
var Rsquare = {
  errtol : 1.0E-010,
  maxiter : 1000
}
  Rsquare.cdf4 = function(paramDouble1, paramDouble2, paramInt, paramDouble3){
    if ((paramInt < 1) || (paramDouble2 < paramInt))  {
      return (0.0 / 0.0);
    }
    if (paramDouble3 < 0.0){
      return (0.0 / 0.0);
    }
    if (paramDouble1 <= 0.0) {
      return 0.0;
    }
    if (paramDouble1 >= 1.0) {
      return 1.0;
    }
    var d1 = paramDouble2 - 1.0;
    var i = (d1 * paramDouble3 / (2.0 * (1.0 - paramDouble3)));
    var d2 = (paramInt - 1.0) / 2.0 + i;
    var d3 = (d1 - paramInt + 1.0) / 2.0;
    var d4 = Beta.cdf3(paramDouble1, d2, d3);
    if (paramDouble3 < 1.0E-012) {
      return d4;
    }
    var d5 = d4;
    var d6 = Math.exp((d2 - 1.0) * Math.log(paramDouble1) + d3 * Math.log(1.0 - paramDouble1) + MoreMath.logGamma(d2 + d3 - 1.0) - MoreMath.logGamma(d2) - MoreMath.logGamma(d3));
    

    var d7 = d6 * (d2 + d3 - 1.0) * paramDouble1 / d2;
    var d8 = Math.exp(MoreMath.logGamma(d1 / 2.0 + i) - MoreMath.logGamma(i + 1) - MoreMath.logGamma(d1 / 2.0) + i * Math.log(paramDouble3) + d1 / 2.0 * Math.log(1.0 - paramDouble3));
    

    var d9 = d8;
    var d10 = 1.0 - d8;
    var d12 = d8 * d4;
    var j = 0;
    for (var k = 1; j == 0; k++){
      d6 *= (d2 + d3 + k - 2.0) * paramDouble1 / (d2 + k - 1.0);
      d4 -= d6;
      d8 = d8 * (d1 / 2.0 + i + k - 1.0) * paramDouble3 / (i + k);
      d12 += d8 * d4;
      var d11 = d10 * d4;
      d10 -= d8;
      if (k > i){
        if ((d11 < this.errtol) || (k > this.maxiter)) {
          j = 1;
        }
      } else {
        d7 *= (d2 - k + 1.0) / (paramDouble1 * (d2 + d3 - k));
        d5 += d7;
        d9 *= (i - k + 1) / (paramDouble3 * (d1 / 2.0 + i - k));
        d12 += d9 * d5;
        d10 -= d9;
        if ((d10 < this.errtol) || (k > this.maxiter)) {
          j = 1;
        }
      }
      if (k > this.maxiter) {
      }
    }
    return d12;
  }
  
Rsquare.cdf3 = function(paramDouble1, paramDouble2, paramInt){
    return this.cdf4(paramDouble1, paramDouble2, paramInt, 0.0);
  }
var RsqAux ={
  xMax : 9.0E+300,
  xMin : -9.0E+300,
  closedMax : false,
  closedMin : false,
  maxIter : 100,
  maxSearch : 25,
  xeps : 1.0E-006,
  feps : 1.0E-010,
  verbose : false,
  N:0.0,
  rho2:0.0,
  p:0
}
RsqAux.news = function(paramDouble1, paramInt, paramDouble2){
  this.N = paramDouble1;
  this.p = paramInt;
  this.rho2 = paramDouble2;
  this.xMin = 0.0;
  this.xMax = 1.0;
  this.closedMin = true;
  this.closedMax = true;
}

RsqAux.of = function(paramDouble){
  return Rsquare.cdf4(paramDouble, this.N, this.p, this.rho2);
}

Rsquare.quantile4 = function (paramDouble1, paramDouble2, paramInt, paramDouble3){
    if ((paramInt < 1) || (paramDouble2 < paramInt))    {
      return (0.0 / 0.0);
    }
    if (paramDouble3 < 0.0)
    {
      return (0.0 / 0.0);
    }
    paramDouble1 = paramDouble1 > 1.0 ? 1.0 : paramDouble1 < 0.0 ? 0.0 : paramDouble1;
    if (paramDouble1 * (1.0 - paramDouble1) == 0.0) {
      return paramInt;
    }
    var d1 = paramDouble3 / (1.0 - paramDouble3);
    var d2 = Beta.quantile(paramDouble1, paramInt, paramDouble2 - paramInt - 1.0, d1 * (paramDouble2 - 1.0));
    RsqAux.news(paramDouble2, paramInt, paramDouble3);
    return Solve.search(RsqAux, paramDouble1, d2, 0.01);
  }
  
Rsquare.quantile3 = function(paramDouble1, paramDouble2, paramInt){
    return this.quantile4(paramDouble1, paramDouble2, paramInt, 0.0);
  }
//end of Rsquare.js

//begin of RsquareGUI.js
var RsquareGUI = {}

RsquareGUI.news = function(){
  this.rho2 = 0.0;
    this.n = 0.0;
    this.alpha = 0.0;
    this.power = 0.0;
    this.preds = 0.0;
}

RsquareGUI.gui = function(){
    this.alpha=0.05;
    this.rho2=0.1;
    this.n=50.0;
    this.preds=1.0;
    this.power=0.0;
}

RsquareGUI.click = function(){
  this.preds = (this.preds < 1.0 ? 1.0 :Math.round(this.preds));
  var i = this.preds + 1;
  this.n =Math.max(this.n, this.preds + 1.0);
  this.rho2 =Math.min(this.rho2, 0.999);
  this.alpha =Math.max(0.0001,Math.min(0.5, this.alpha));
  var d = Rsquare.quantile3(1.0 - this.alpha, this.n, i);
  this.power = (1.0 - Rsquare.cdf4(d, this.n, i, this.rho2));
}

module.exports.RsquareGUI_handle = function(id, key, value){
     if(isNullOrEmpty(RsquareGUI.n)){
      RsquareGUI.news(); 
         RsquareGUI.gui();
         RsquareGUI.click();
     }
     if(!isNullOrEmpty(id)){
          if(!isNullOrEmpty(key)){
            if(key=="alpha"){
              if(!isNullOrEmpty(value)){
                RsquareGUI.alpha = parseFloat(value);
                RsquareGUI.click();
              }
            }
            if(key=="rho2"){
              if(!isNullOrEmpty(value)){
                RsquareGUI.rho2 = parseFloat(value);
                RsquareGUI.click();
              }
            }
            if(key=="n"){
              if(!isNullOrEmpty(value)){
                RsquareGUI.n = parseFloat(value);
                RsquareGUI.click();
              }
            }
            if(key=="preds"){
              if(!isNullOrEmpty(value)){
                RsquareGUI.preds = parseFloat(value);
                RsquareGUI.click();
                }
            }
            if(key=="power"){
              if(!isNullOrEmpty(value)){
                RsquareGUI.power = parseFloat(value);
                RsquareGUI.click();
              }
            }
          }
     }else{
//       RsquareGUI.news();
//       RsquareGUI.gui();  
     }
     return RsquareGUI;
}
//end of RsquareGUI.js

//Begin of SimpleChi2GUI.js
var SimpleChi2GUI={}
SimpleChi2GUI.news = function(){
  this.proChi2 = 0.0;
  this.proN = 0.0;
  this.n = 0.0;
  this.df = 0.0;
  this.Alpha = 0.0;
  this.Power = 0.0;
}
SimpleChi2GUI.gui = function() {
    this.proChi2=10.0;
    this.proN=100.0;
    this.df=3.0;
    this.Alpha=0.05;
    this.n=50.0;
    this.Power=0.6;
}
    
SimpleChi2GUI.click = function ()   {
    this.n = Math.max(Math.round(this.n), 2.0);
    this.df = Math.max(Math.round(this.df), 1.0);
    this.Power = Chi2.power3(this.n * this.proChi2 / this.proN, this.df, this.Alpha);
  }
    
SimpleChi2GUI.Power_changed = function ()   {
    this.n = (Chi2.lambda(this.Power, this.df, this.Alpha) * this.proN / this.proChi2);
    this.click();
  }

module.exports.SimpleChi2GUI_handle = function(id, key, value){
   if(isNullOrEmpty(SimpleChi2GUI.n)){
      SimpleChi2GUI.news(); 
         SimpleChi2GUI.gui();
         SimpleChi2GUI.click();
     }
     if(!isNullOrEmpty(id)){
          if(!isNullOrEmpty(key)){
            if(key=="proChi2"){
              if(!isNullOrEmpty(value)){
                SimpleChi2GUI.proChi2 = parseFloat(value);
                SimpleChi2GUI.click();
              }
            }
            if(key=="proN"){
              if(!isNullOrEmpty(value)){
                 SimpleChi2GUI.proN = parseFloat(value);
                 SimpleChi2GUI.click();
              }
             }
             if(key=="n"){
               if(!isNullOrEmpty(value)){
                 SimpleChi2GUI.n = parseFloat(value);
                 SimpleChi2GUI.click();
               }
             }
             if(key=="df"){
               if(!isNullOrEmpty(value)){
                 SimpleChi2GUI.df = parseFloat(value);
                 SimpleChi2GUI.click();
               }
             }
             if(key=="Alpha"){
               if(!isNullOrEmpty(value)){
                 SimpleChi2GUI.Alpha = parseFloat(value);
                 SimpleChi2GUI.click();
               }
             }
             if(key=="Power"){
               if(!isNullOrEmpty(value)){
                 SimpleChi2GUI.Power = parseFloat(value);
                 SimpleChi2GUI.Power_changed();
               }
             }
          }
     }else{
       SimpleChi2GUI.news();
       SimpleChi2GUI.gui(); 
     }
     return SimpleChi2GUI;
}
//end of SimpleChi2GUI

//begin of SimplePoissonGUI.js
var SimplePoissonGUI = {}

SimplePoissonGUI.news = function(){
  this.lambda0 = 0.0;
  this.lambda = 0.0;
  this.alpha = 0.0;
  this.power = 0.0;
  this.sizes = 0.0;
  this.lower = 0.0;
  this.upper = 0.0;
  this.n = 0.0;
  this.alt = 0;
}

SimplePoissonGUI.gui = function (){
  this.lambda0=1.0;
    this.alt=1;
    this.alpha=0.05;
    this.lower=0.0;
    this.upper=0.0;
    this.sizes=0.0;
    this.lambda=1.0;
    this.n=50.0;
    this.power=0.0;
}

SimplePoissonGUI.click = function (){
  var d = this.alt == 1 ? this.alpha / 2.0 : this.alpha;
  this.sizes = (this.power = 0.0);
  if (this.alt < 2)  {
    this.lower = (Poisson.quantile(d, this.n * this.lambda0) + 1);
    this.sizes += Poisson.cdf(this.lower - 1, this.n * this.lambda0);
    this.power += Poisson.cdf(this.lower - 1, this.n * this.lambda);
  }
  if (this.alt > 0)  {
    this.upper = Poisson.quantile(1.0 - d, this.n * this.lambda0);
    this.sizes += 1.0 - Poisson.cdf(this.upper + 1, this.n * this.lambda0);
    this.power += 1.0 - Poisson.cdf(this.upper + 1, this.n * this.lambda);
  }
}

SimplePoissonGUI.alt_changed = function (){
    setVisible("lower", this.alt < 2);
    setVisible("upper", this.alt > 0);
    this.click();
}

module.exports.SimplePoissonGUI_handle = function(id, key, value){
  if(isNullOrEmpty(SimplePoissonGUI.n)){
    SimplePoissonGUI.news(); 
    SimplePoissonGUI.gui();
    SimplePoissonGUI.click();
    }
    if(!isNullOrEmpty(id)){
          if(!isNullOrEmpty(key)){
            if(key=="lambda0"){
              if(!isNullOrEmpty(value)){
                SimplePoissonGUI.lambda0 = parseFloat(value);
                SimplePoissonGUI.click();
              }
            }
            if(key=="lambda"){
              if(!isNullOrEmpty(value)){
                SimplePoissonGUI.lambda = parseFloat(value);
                SimplePoissonGUI.click();
              }
            }
            if(key=="alpha"){
              if(!isNullOrEmpty(value)){
                SimplePoissonGUI.alpha = parseFloat(value);
                SimplePoissonGUI.click();
              }
            }
            if(key=="n"){
                if(!isNullOrEmpty(value)){
                  SimplePoissonGUI.n = parseFloat(value);
                  SimplePoissonGUI.click();
                }
            }
            if(key=="alt"){
              if(!isNullOrEmpty(value)){
                  SimplePoissonGUI.alt = parseInt(value);
                SimplePoissonGUI.alt_changed();
              }
            }
          }
    }else{
      SimplePoissonGUI.news();
      SimplePoissonGUI.gui(); 
      SimplePoissonGUI.click();
    }
    return SimplePoissonGUI;
}
//end of SimplePoissonGUI.js

//begin of TwoTGUI.js
function TwoTGUI(){}

TwoTGUI.prototype = new PowerCalculator();
TwoTGUI.prototype.sigma1 = 0.0;
TwoTGUI.prototype.sigma2 = 0.0;
TwoTGUI.prototype.n1 = 0.0;
TwoTGUI.prototype.n2 = 0.0;
TwoTGUI.prototype.diff = 0.0;
TwoTGUI.prototype.thresh = 0.0;
TwoTGUI.prototype.alpha = 0.0;
TwoTGUI.prototype.df = 0.0;
TwoTGUI.prototype.power = 0.0;
TwoTGUI.prototype.v1 = 0.0;
TwoTGUI.prototype.v2 = 0.0;
TwoTGUI.prototype.delta = 0.0;
TwoTGUI.prototype.mult = 0.0;
TwoTGUI.prototype.saveN1 = 0.0;
TwoTGUI.prototype.saveN2 = 0.0;
TwoTGUI.prototype.eqs = 0;
TwoTGUI.prototype.alloc = 0;
TwoTGUI.prototype.tt = 0;
TwoTGUI.prototype.prevTT = 0;
TwoTGUI.prototype.opt = 0;
TwoTGUI.prototype.equiv = 0;
TwoTGUI.prototype.rocMeth = 0;

TwoTGUI.prototype.news = function(){
    this.sigma1 = 0.0;
    this.sigma2 = 0.0;
    this.n1 = 0.0;
    this.n2 = 0.0;
    this.diff = 0.0;
    this.thresh = 0.0;
    this.alpha = 0.0;
    this.df = 0.0;
    this.power = 0.0;
    this.v1 = 0.0;
    this.v2 = 0.0;
    this.delta = 0.0;
    this.mult = 0.0;
    this.saveN1 = 0.0;
    this.saveN2 = 0.0;
    this.eqs = 0;
    this.alloc = 0;
    this.tt = 0;
    this.prevTT = 0;
    this.opt = 0;
    this.equiv = 0;
    this.rocMeth = 0;
  }

TwoTGUI.prototype.gui = function(){
  this.sigma1=1.0;
  this.sigma2=1.0;
  this.eqs=1;
  this.n1=25.0;
  this.n2=25.0;
  this.alloc=1;
  this.tt=1;
    this.alpha=0.05;
  this.equiv=0;
  this.thresh=1.0;
  this.df=this.n1 + this.n2 - 2.0;
  this.diff=0.5;
  this.power=0.5;
  this.opt=0;
  this.prevTT = this.tt;
}

TwoTGUI.prototype.click = function(){
  this.n1 = Math.max(2.0,Math.round(this.n1));
  if (this.equiv == 1) {
    this.tt = 1;
  }
  if (this.eqs == 1) {
    this.sigma2 = this.sigma1;
  }
  if (this.alloc == 1) {
    this.n2 = this.n1;
  } else if (this.alloc == 2) {
    this.n2 =Math.max(2.0,Math.round(this.n1 * this.sigma2 / this.sigma1));
  }
  this.sattPower();
}

TwoTGUI.prototype.sigma2_changed = function(){
  if (this.eqs == 1) {
    this.sigma1 = this.sigma2;
  }
  if (this.alloc == 2) {
    this.n2 =Math.max(2.0,Math.round(this.n1 * this.sigma2 / this.sigma1));
  }
  this.sattPower();
}

TwoTGUI.prototype.n2_changed = function(){
  this.n2 =Math.max(2.0,Math.round(this.n2));
  if (this.alloc == 1) {
    this.n1 = this.n2;
  } else if (this.alloc == 2) {
    this.n1 =Math.max(2.0,Math.round(this.n2 * this.sigma1 / this.sigma2));
  }
  this.sattPower();
}

TwoTGUI.prototype.sattPower = function(){
  this.v1 = (this.sigma1 * this.sigma1 / this.n1);
  this.v2 = (this.sigma2 * this.sigma2 / this.n2);
  var d =Math.sqrt(this.v1 + this.v2);
  this.delta = (this.diff / d);
  this.df = (this.eqs == 1 ? this.n1 + this.n2 - 2.0 : (this.v1 + this.v2) * (this.v1 + this.v2) / (this.v1 * this.v1 / (this.n1 - 1.0) + this.v2 * this.v2 / (this.n2 - 1.0)));
  if (this.equiv == 0) {
    this.power = (this.rocMeth == 0 ? T.power(this.delta, this.df, 1 - this.tt, this.alpha) : T.rocArea(this.delta, this.df, 1 - this.tt));
  } else {
    this.power = (this.rocMeth == 0 ? T.powerEquiv5(this.diff, this.thresh, d, this.df, this.alpha) : T.rocEquiv(this.diff, this.thresh, d, this.df));
  }
}

TwoTGUI.prototype.power_changed = function(){
  if ((this.equiv == 1) || (this.rocMeth == 1)){
    this.power_changed_numerical();
    return;
  }
  switch (this.opt){
  case 0: 
    var d = this.power;
    this.diff =Math.max(this.diff, 0.01 * (this.sigma1 + this.sigma2));
    for (var i = 0; i < 3; i++){
      this.delta = T.delta(d, this.df, 1 - this.tt, this.alpha);
      this.mult = ((this.v1 + this.v2) * this.delta * this.delta / this.diff / this.diff);
      this.n1 *= this.mult;
      this.n2 *= this.mult;
      this.sattPower();
    }
    this.n1 =Math.max(Math.round(this.n1), 2.0);
    this.n2 =Math.max(Math.round(this.n2), 2.0);
    this.sattPower();
    break;
  case 1: 
    this.delta = T.delta(this.power, this.df, 1 - this.tt, this.alpha);
    this.diff = (this.delta *Math.sqrt(this.v1 + this.v2));
    this.sattPower();
  }
}

TwoTGUI.prototype.power_changed_numerical = function(){
  switch (this.opt){
  case 0: 
    this.saveN1 = this.n1;
    this.saveN2 = this.n2;
    PowerCalculatorAux.news("mult", "power", this);
    PowerCalculatorAux.xMin =Math.max(2.0 / this.n1, 2.0 / this.n2);
    PowerCalculatorAux.closedMin = true;
    PowerCalculatorAux.xeps =Math.max(0.5 / this.n1, 0.5 / this.n2);
    this.mult =this.solve(PowerCalculatorAux, this.power, 1.0, 0.1);
    this.n1 = Math.max(Math.round(this.n1 * this.mult), 10.0);
    this.n2 = Math.max(Math.round(this.n2 * this.mult), 10.0);
    break;
  case 1: 
    PowerCalculatorAux.news("diff", "power", this);
    PowerCalculatorAux.xeps = (0.005 * (this.sigma1 + this.sigma2));
    this.diff = Math.max(this.diff, 0.1 * (this.sigma1 + this.sigma2));
    this.diff =this.solve(localPowerCalculatorAux2, this.power, this.diff, 0.1 * this.diff);
  }
  this.sattPower();
}

TwoTGUI.prototype.mult_changed = function(){
  this.n1 = (this.saveN1 * this.mult);
  this.n2 = (this.saveN2 * this.mult);
  this.sattPower();
}



TwoTGUI.prototype.equiv_changed = function(){
  if (this.equiv == 1){
    this.prevTT = this.tt;
    this.tt = 1;
  }else{
    this.tt = this.prevTT;
  }
  this.sattPower();
}

TwoTGUI.prototype.rocMeth_changed = function(){
  if (this.rocMeth == 1) {
  } else {
  }
 //setVisible("alpha", this.rocMeth == 0);
  this.sattPower();
}

var ttg = null;
module.exports.TwoTGUI_dataDrivenMode_getpower = function(sig1_in, sig2_in, n1_in, n2_in, alpha_in, df_in, diff_in){
  ttg = new TwoTGUI();
  ttg.sigma1 = parseFloat(sig1_in);
  ttg.sigma2 = parseFloat(sig2_in);
  ttg.n1 = parseFloat(n1_in);
  ttg.n2 = parseFloat(n2_in);
  ttg.alpha = parseFloat(alpha_in);
  ttg.diff = parseFloat(diff_in);
  ttg.eqs = 0;
  ttg.alloc = 0;
  ttg.tt = 0;
  ttg.equiv = 0;
  ttg.click();
  return ttg;
}
module.exports.TwoTGUI_dataDrivenMode_changepower = function(sig1_in, sig2_in, n1_in, n2_in, alpha_in, df_in, diff_in, power_in){
  ttg = new TwoTGUI();
  ttg.sigma1 = parseFloat(sig1_in);
  ttg.sigma2 = parseFloat(sig2_in);
  ttg.n1 = parseFloat(n1_in);
  ttg.n2 = parseFloat(n2_in);
  ttg.alpha = parseFloat(alpha_in);
  ttg.diff = parseFloat(diff_in);
  ttg.eqs = 0;
  ttg.alloc = 0;
  ttg.tt = 0;
  ttg.equiv = 0;
  ttg.click();
  ttg.power = parseFloat(power_in)
  var temp_power = parseFloat(power_in)
  ttg.power_changed()
  if(!ttg.n1){
    ttg.n1 = 2;
    ttg.n2 = ttg.n1;
  }
  if(!ttg.power){
    ttg.power = temp_power;
  }
  if(!ttg.df){
    ttg.df = 1;
  }
  return ttg;
}



module.exports.TwoTGUI_handle = function(id, key, value){
      if(isNullOrEmpty(ttg)){
        ttg = new TwoTGUI();
        ttg.news(); 
          ttg.gui();
          ttg.click();
      }
      if(!isNullOrEmpty(id)){
            if(key=="sigma1"){
              if(!isNullOrEmpty(value)){
                ttg.sigma1 = parseFloat(value);
                ttg.click();
              }
            }
            if(key=="sigma2"){
              if(!isNullOrEmpty(value)){
                ttg.sigma2 = parseFloat(value);
                ttg.sigma2_changed();
              }
            }
            if(key=="eqs"){
               if(isNullOrEmpty(value)){
                 ttg.eqs = 0;
               }else{
                 ttg.eqs = 1;
               }
               ttg.click();
            }
            if(key=="n1"){
              if(!isNullOrEmpty(value)){
                ttg.n1 = parseFloat(value);
                ttg.click();
              }
            }
            if(key=="n2"){
              if(!isNullOrEmpty(value)){
                ttg.n2 = parseFloat(value);
                ttg.n2_changed();
              }
            }
            if(key=="alloc"){
              if(!isNullOrEmpty(value)){
                ttg.alloc = parseInt(value);
                ttg.click();
              }
            }
            if(key=="tt"){
              if(isNullOrEmpty(value)){
                ttg.tt = 0;
              }else{
                ttg.tt = 1;
              }
              ttg.click();
            }
            if(key=="alpha"){
              if(!isNullOrEmpty(value)){
                ttg.alpha = parseFloat(value);
                ttg.click();
              }
            }
            if(key=="equiv"){
              if(isNullOrEmpty(value)){
                ttg.equiv = 0;
              }else{
                ttg.equiv = 1;
              }
              ttg.equiv_changed();
            }
            if(key=="thresh"){
               if(!isNullOrEmpty(value)){
                 ttg.thresh = parseFloat(value);
                 ttg.click();
               }
            }
            if(key=="diff"){
              if(!isNullOrEmpty(value)){
                ttg.diff = parseFloat(value);
                ttg.click();
              }
            }
            if(key=="power"){
               if(!isNullOrEmpty(value)){
                 ttg.power = parseFloat(value);
                 ttg.power_changed();
               }
            }
            if(key=="opt"){
                if(!isNullOrEmpty(value)){
                  ttg.opt = parseInt(value);
                  ttg.click();
                }
            }
      }else{
        ttg.news();
        ttg.gui();  
      }
      return ttg;
}
//end of TwoTGUI.js
% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:36
% DurationCPUTime: 5.06s
% Computational Cost: add. (9017->330), mult. (10069->428), div. (0->0), fcn. (10246->8), ass. (0->184)
t248 = qJ(1) + pkin(8);
t232 = cos(t248);
t145 = t232 * qJ(3);
t231 = sin(t248);
t226 = t231 * pkin(2);
t112 = t226 - t145;
t142 = qJD(3) * t231;
t211 = qJD(1) * t232;
t256 = qJ(3) * t211 + t142;
t327 = -qJD(1) * t112 + t142 - t256;
t156 = cos(qJ(5));
t154 = sin(qJ(5));
t294 = sin(qJ(4));
t295 = cos(qJ(4));
t108 = -t231 * t294 - t232 * t295;
t109 = -t231 * t295 + t232 * t294;
t274 = Icges(6,4) * t156;
t207 = -Icges(6,2) * t154 + t274;
t60 = -Icges(6,6) * t109 + t108 * t207;
t279 = t154 * t60;
t275 = Icges(6,4) * t154;
t209 = Icges(6,1) * t156 - t275;
t64 = -Icges(6,5) * t109 + t108 * t209;
t212 = t156 * t64 - t279;
t61 = Icges(6,6) * t108 + t109 * t207;
t280 = t154 * t61;
t65 = Icges(6,5) * t108 + t109 * t209;
t326 = -t156 * t65 + t280;
t262 = t109 * t156;
t205 = Icges(6,5) * t156 - Icges(6,6) * t154;
t57 = Icges(6,3) * t108 + t109 * t205;
t289 = -t108 * t57 - t262 * t65;
t56 = -Icges(6,3) * t109 + t108 * t205;
t325 = -(t56 - t280) * t109 + t289;
t252 = qJD(5) * t109;
t253 = qJD(5) * t108;
t153 = qJD(1) - qJD(4);
t296 = -t153 / 0.2e1;
t36 = -t154 * t65 - t156 * t61;
t37 = -t154 * t64 - t156 * t60;
t324 = ((-t37 * t108 + t36 * t109) * qJD(5) - t36 * t252 + t37 * t253) * t296;
t322 = t109 * t56;
t317 = t57 * t109;
t315 = t57 + t279;
t150 = t232 * pkin(3);
t152 = cos(qJ(1)) * pkin(1);
t254 = t232 * pkin(2) + t231 * qJ(3);
t163 = -t254 - t152;
t229 = t150 - t163;
t175 = t232 * rSges(4,1) + t231 * rSges(4,3);
t306 = t163 - t175;
t143 = qJD(3) * t232;
t314 = -qJD(1) * t229 + t143;
t291 = sin(qJ(1)) * pkin(1);
t185 = -pkin(3) * t231 - t291;
t166 = -t226 + t185;
t312 = -t327 + (t166 - t185) * qJD(1);
t220 = rSges(6,1) * t156 - rSges(6,2) * t154;
t119 = t220 * qJD(5);
t158 = qJD(1) ^ 2;
t210 = qJD(1) * t231;
t276 = (-pkin(2) * t210 + t142 + t256) * qJD(1);
t173 = t158 * t185 + t276;
t219 = rSges(6,1) * t154 + rSges(6,2) * t156;
t93 = t153 * t108;
t94 = t153 * t109;
t243 = t94 * pkin(4) + t93 * pkin(7);
t250 = qJD(5) * t154;
t195 = t108 * t250 + t156 * t94;
t192 = t195 * rSges(6,1) + t93 * rSges(6,3);
t249 = qJD(5) * t156;
t196 = t108 * t249 - t154 * t94;
t30 = rSges(6,2) * t196 + t192;
t11 = (t243 + t30) * t153 + (t109 * t119 + t219 * t93) * qJD(5) + t173;
t186 = -t150 - t152;
t230 = (-qJD(1) * t254 + 0.2e1 * t143) * qJD(1);
t162 = t158 * t186 + t230;
t244 = t93 * pkin(4) - t94 * pkin(7);
t193 = t109 * t250 - t156 * t93;
t191 = t193 * rSges(6,1) + t94 * rSges(6,3);
t194 = t109 * t249 + t154 * t93;
t29 = rSges(6,2) * t194 + t191;
t12 = (t244 - t29) * t153 + (t108 * t119 - t219 * t94) * qJD(5) + t162;
t172 = t142 + (-t112 + t185) * qJD(1);
t251 = qJD(5) * t219;
t105 = t108 * rSges(6,3);
t260 = rSges(6,1) * t262 + t105;
t263 = t109 * t154;
t70 = rSges(6,2) * t263 - t260;
t97 = -t109 * pkin(4) - t108 * pkin(7);
t303 = t108 * t251 - (t70 + t97) * t153;
t27 = t172 + t303;
t171 = (t254 - t186) * qJD(1) - t143;
t233 = -t108 * pkin(4) + pkin(7) * t109;
t240 = t109 * t251;
t104 = t109 * rSges(6,3);
t265 = t108 * t156;
t259 = -rSges(6,1) * t265 + t104;
t266 = t108 * t154;
t71 = rSges(6,2) * t266 + t259;
t28 = t240 + (t233 + t71) * t153 + t171;
t311 = (t154 * (t108 * t11 - t109 * t12 - t27 * t93 - t28 * t94) + (t108 * t28 - t109 * t27) * t249) * rSges(6,2);
t222 = rSges(5,1) * t109 - rSges(5,2) * t108;
t310 = t153 * t222;
t204 = Icges(6,5) * t154 + Icges(6,6) * t156;
t80 = t204 * t108;
t79 = t204 * t109;
t309 = 0.2e1 * qJD(5);
t206 = Icges(6,2) * t156 + t275;
t208 = Icges(6,1) * t154 + t274;
t203 = -t154 * t206 + t156 * t208;
t40 = t109 * t203 + t80;
t38 = t40 * t153;
t41 = t108 * t203 - t79;
t39 = t41 * t153;
t308 = rSges(3,1) * t232 - rSges(3,2) * t231;
t305 = -rSges(4,1) * t231 - t291;
t283 = t109 * t206 - t65;
t285 = -t109 * t208 - t61;
t302 = t154 * t283 + t156 * t285;
t288 = t108 * t56 + t262 * t64;
t287 = t265 * t65 - t317;
t286 = t265 * t64 - t322;
t284 = -t108 * t208 - t60;
t282 = t108 * t206 - t64;
t261 = t205 * t153;
t258 = -t206 + t209;
t257 = -t207 - t208;
t247 = t158 * t291;
t246 = t158 * t152;
t236 = t253 / 0.2e1;
t235 = -t252 / 0.2e1;
t228 = -t232 / 0.2e1;
t227 = t231 / 0.2e1;
t49 = rSges(5,1) * t94 - rSges(5,2) * t93;
t48 = -rSges(5,1) * t93 - rSges(5,2) * t94;
t221 = rSges(5,1) * t108 + rSges(5,2) * t109;
t217 = -t108 * t27 - t109 * t28;
t216 = t108 * t71 + t109 * t70;
t202 = -t97 + t260;
t201 = -t233 - t259;
t17 = -t266 * t61 + t287;
t18 = -t266 * t60 + t286;
t189 = qJD(5) * (-t108 * t17 + t109 * t18);
t15 = -t263 * t61 - t289;
t16 = -t263 * t60 + t288;
t187 = (-t108 * t15 + t109 * t16) * qJD(5);
t184 = t154 * t282 + t156 * t284;
t183 = -t192 - t243;
t182 = t191 - t244;
t180 = t108 * t30 + t109 * t29 + t70 * t93 - t71 * t94;
t179 = (t154 * t257 + t156 * t258) * t153;
t114 = rSges(3,1) * t231 + rSges(3,2) * t232;
t117 = t207 * qJD(5);
t118 = t209 * qJD(5);
t167 = -t117 * t154 + t118 * t156 + (-t154 * t208 - t156 * t206) * qJD(5);
t165 = -t226 + t305;
t161 = t145 + t166;
t10 = qJD(5) * t212 - t154 * (Icges(6,1) * t195 + Icges(6,4) * t196 + Icges(6,5) * t93) - t156 * (Icges(6,4) * t195 + Icges(6,2) * t196 + Icges(6,6) * t93);
t116 = t205 * qJD(5);
t13 = t108 * t116 + t109 * t167 + t203 * t93 - t204 * t94;
t14 = t108 * t167 - t109 * t116 - t203 * t94 - t204 * t93;
t5 = t38 + t187;
t6 = t39 + t189;
t9 = -qJD(5) * t326 - t154 * (Icges(6,1) * t193 + Icges(6,4) * t194 + Icges(6,5) * t94) - t156 * (Icges(6,4) * t193 + Icges(6,2) * t194 + Icges(6,6) * t94);
t159 = t153 * (-t118 * t154 + t206 * t250 + (-qJD(5) * t208 - t117) * t156) + t5 * t252 / 0.2e1 + (t14 + t10) * t235 + (t13 + t6 + t9) * t236 - ((t41 - t37) * t93 + (t40 - t36) * t94) * qJD(5) / 0.2e1;
t147 = t232 * rSges(4,3);
t141 = rSges(4,3) * t211;
t86 = t219 * t108;
t85 = t219 * t109;
t74 = -t158 * t175 + t230 - t246;
t73 = -t247 + qJD(1) * (-rSges(4,1) * t210 + t141) + t276;
t69 = t109 * t220 + t105;
t68 = t108 * t220 - t104;
t55 = -t153 * t221 + t171;
t54 = t172 + t310;
t35 = -t153 * t48 + t162;
t34 = t153 * t49 + t173;
t31 = qJD(5) * t216 + qJD(2);
t22 = Icges(6,5) * t195 + Icges(6,6) * t196 + Icges(6,3) * t93;
t21 = Icges(6,5) * t193 + Icges(6,6) * t194 + Icges(6,3) * t94;
t8 = (-t108 * t61 - t109 * t60) * t154 + t287 + t288;
t7 = t180 * qJD(5);
t1 = [-t159 + t324 + m(3) * ((-t158 * t308 - t246) * (-t114 - t291) + (-t114 * t158 - t247) * (t152 + t308)) + (-t38 + ((t17 + (-t212 + t57) * t109) * t109 + (t315 * t108 + t18 - t286 - t322) * t108) * qJD(5)) * t235 + (t39 + ((t109 * t326 + t15 + t286) * t109 + (-t109 * t315 + t16 - t8) * t108) * qJD(5)) * t236 + (t12 * (t161 + t202) + t11 * (-t201 + t229) + t311 - t31 * (t69 + t70) * t253 + (-t182 + t314) * t27 + (-t183 + t27 - (t69 - t97) * t153 - t219 * t253 + t312) * t28) * m(6) + (t35 * (t161 + t222) + t34 * (-t221 + t229) + (t314 - t48) * t54 + (-t310 + t312 + t49 + t54) * t55) * m(5) + (t74 * (t145 + t147 + t165) - t73 * t306 + (-t141 - (t165 - t147 - t305) * qJD(1) + t327) * (qJD(1) * t306 + t143)) * m(4); m(6) * t7; 0.2e1 * (t11 * t228 + t12 * t227) * m(6) + 0.2e1 * (t227 * t35 + t228 * t34) * m(5) + 0.2e1 * (t227 * t74 + t228 * t73) * m(4); t159 + t324 + (t38 + ((-t17 + t8) * t109 + (t108 * t212 - t18 + t325) * t108) * qJD(5)) * t235 + (-t39 + ((-t15 - t325) * t109 + (-t16 + (-t326 + t56) * t108 - t317) * t108) * qJD(5)) * t236 + (-t12 * t202 + t11 * t201 - t311 - t31 * (t68 + t71) * t252 + (t183 + t303) * t28 + (t182 - t240 - (-t68 + t233) * t153) * t27) * m(6) + (-(-t221 * t54 - t222 * t55) * t153 + t34 * t221 - t35 * t222 + t48 * t54 - t49 * t55) * m(5); t153 * (t10 * t109 - t108 * t9 - t36 * t94 - t37 * t93) / 0.2e1 + ((t80 * t252 - t261) * t109 + (t179 + (-t302 * t108 + (-t79 + t184) * t109) * qJD(5)) * t108) * t235 + ((t79 * t253 + t261) * t108 + (t179 + (t184 * t109 + (-t302 - t80) * t108) * qJD(5)) * t109) * t236 + ((t154 * t258 - t156 * t257) * t153 + ((t108 * t283 - t109 * t282) * t156 + (-t108 * t285 + t109 * t284) * t154) * qJD(5)) * t296 - (t13 * t153 + (-(-t108 * t21 - t326 * t93 - t57 * t94) * t108 + t109 * (-t108 * t22 + t212 * t93 - t56 * t94) + t15 * t94 + t16 * t93) * t309) * t108 / 0.2e1 + (t14 * t153 + (-t108 * (t109 * t21 + t326 * t94 - t57 * t93) + t109 * (t109 * t22 - t212 * t94 - t56 * t93) + t17 * t94 + t18 * t93) * t309) * t109 / 0.2e1 + (t6 + t189) * t93 / 0.2e1 + (t5 + t187) * t94 / 0.2e1 + (t7 * t216 + t31 * t180 - t217 * t119 - (-t108 * t12 - t109 * t11 + t27 * t94 - t28 * t93) * t219 - (-t27 * t85 + t28 * t86) * t153 - (t31 * (t108 * t86 + t109 * t85) - t217 * t220) * qJD(5)) * m(6);];
tauc = t1(:);

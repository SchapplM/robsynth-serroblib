% Calculate time derivative of joint inertia matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:21
% EndTime: 2019-12-31 18:12:33
% DurationCPUTime: 7.59s
% Computational Cost: add. (3970->407), mult. (5400->575), div. (0->0), fcn. (4114->6), ass. (0->207)
t153 = pkin(7) + qJ(3);
t145 = sin(t153);
t146 = cos(t153);
t263 = Icges(5,6) * t146;
t271 = Icges(4,4) * t146;
t310 = t263 + t271 + (Icges(4,1) + Icges(5,2)) * t145;
t245 = qJD(3) * t146;
t262 = Icges(6,6) * t145;
t264 = Icges(5,6) * t145;
t272 = Icges(4,4) * t145;
t309 = -t262 + t264 + t272 + (Icges(4,2) + Icges(6,2) + Icges(5,3)) * t146;
t278 = rSges(6,3) + qJ(5);
t308 = t145 * t278;
t159 = sin(qJ(1));
t288 = rSges(6,1) + pkin(4);
t299 = t159 * t288;
t307 = t310 * qJD(3);
t306 = t309 * t245;
t305 = t159 / 0.2e1;
t160 = cos(qJ(1));
t304 = -t160 / 0.2e1;
t303 = -qJD(1) / 0.2e1;
t302 = qJD(1) / 0.2e1;
t243 = qJD(3) * t160;
t231 = t145 * t243;
t248 = qJD(1) * t159;
t301 = t146 * t248 + t231;
t193 = Icges(5,2) * t146 - t264;
t295 = Icges(5,4) * t160 + t159 * t193;
t189 = -Icges(5,3) * t145 + t263;
t296 = Icges(5,5) * t160 + t159 * t189;
t205 = -t145 * t296 + t146 * t295;
t177 = t205 * t160;
t70 = Icges(5,5) * t159 - t160 * t189;
t74 = Icges(5,4) * t159 - t160 * t193;
t206 = t145 * t70 - t146 * t74;
t178 = t206 * t159;
t186 = Icges(6,3) * t146 + t262;
t69 = -Icges(6,5) * t160 + t159 * t186;
t261 = Icges(6,6) * t146;
t190 = Icges(6,2) * t145 + t261;
t73 = -Icges(6,4) * t160 + t159 * t190;
t207 = t145 * t73 + t146 * t69;
t179 = t207 * t160;
t68 = Icges(6,5) * t159 + t160 * t186;
t72 = Icges(6,4) * t159 + t160 * t190;
t208 = t145 * t72 + t146 * t68;
t180 = t208 * t159;
t198 = -Icges(4,2) * t145 + t271;
t65 = Icges(4,6) * t159 + t160 * t198;
t200 = Icges(4,1) * t146 - t272;
t67 = Icges(4,5) * t159 + t160 * t200;
t209 = t145 * t65 - t146 * t67;
t181 = t209 * t159;
t64 = -Icges(4,6) * t160 + t159 * t198;
t66 = -Icges(4,5) * t160 + t159 * t200;
t210 = t145 * t64 - t146 * t66;
t182 = t210 * t160;
t151 = t159 * rSges(5,1);
t254 = t146 * t160;
t298 = -rSges(5,2) * t254 + t151;
t149 = t159 * rSges(4,3);
t256 = t145 * t160;
t297 = -rSges(4,2) * t256 + t149;
t194 = Icges(4,5) * t146 - Icges(4,6) * t145;
t62 = -Icges(4,3) * t160 + t159 * t194;
t195 = Icges(6,4) * t145 + Icges(6,5) * t146;
t77 = -Icges(6,1) * t160 + t159 * t195;
t196 = Icges(5,4) * t146 - Icges(5,5) * t145;
t294 = Icges(5,1) * t160 + t159 * t196;
t157 = cos(pkin(7));
t185 = rSges(3,1) * t157 - rSges(3,2) * sin(pkin(7)) + pkin(1);
t280 = rSges(3,3) + qJ(2);
t61 = t159 * t280 + t160 * t185;
t293 = 2 * m(4);
t292 = 2 * m(5);
t291 = 2 * m(6);
t290 = m(5) / 0.2e1;
t289 = m(6) / 0.2e1;
t118 = rSges(4,1) * t145 + rSges(4,2) * t146;
t287 = m(4) * t118;
t286 = pkin(3) * t145;
t285 = pkin(3) * t146;
t217 = qJ(4) * t145 + t285;
t90 = t217 * t159;
t136 = pkin(3) * t254;
t91 = qJ(4) * t256 + t136;
t284 = t159 * t90 + t160 * t91;
t283 = rSges(4,1) * t146;
t282 = rSges(5,2) * t145;
t281 = -rSges(6,2) - qJ(4);
t279 = -rSges(5,3) - qJ(4);
t220 = -rSges(5,2) * t146 + rSges(5,3) * t145;
t81 = qJD(3) * t217 - qJD(4) * t146;
t277 = -t220 * qJD(3) - t81;
t276 = rSges(6,2) * t256 + t254 * t278 + t299;
t218 = rSges(6,2) * t145 + rSges(6,3) * t146;
t255 = t146 * t159;
t275 = qJ(5) * t255 + t159 * t218 - t160 * t288;
t63 = Icges(4,3) * t159 + t160 * t194;
t259 = qJD(1) * t63;
t76 = Icges(6,1) * t159 + t160 * t195;
t258 = qJD(1) * t76;
t78 = Icges(5,1) * t159 - t160 * t196;
t257 = qJD(1) * t78;
t158 = -pkin(6) - qJ(2);
t253 = t158 * t160;
t116 = -qJ(4) * t146 + t286;
t219 = rSges(5,3) * t146 + t282;
t252 = -t116 + t219;
t230 = t146 * t243;
t242 = qJD(4) * t145;
t251 = qJ(4) * t230 + t160 * t242;
t148 = qJD(2) * t160;
t250 = t158 * t248 + t148;
t249 = t159 ^ 2 + t160 ^ 2;
t247 = qJD(1) * t160;
t246 = qJD(3) * t145;
t244 = qJD(3) * t159;
t241 = qJD(5) * t146;
t240 = -pkin(3) - t278;
t239 = -0.1e1 + t249;
t126 = t244 * t286;
t233 = t145 * t248;
t238 = t159 * (qJD(1) * t136 + t159 * t242 - t126 + (t145 * t247 + t146 * t244) * qJ(4)) + t160 * (-pkin(3) * t301 - qJ(4) * t233 + t251) + t90 * t247;
t147 = qJD(2) * t159;
t235 = t147 + t251;
t234 = t126 + t250;
t59 = t252 * t160;
t229 = -rSges(6,2) * t146 + t308;
t138 = pkin(2) * t157 + pkin(1);
t228 = t160 * t138 - t159 * t158;
t227 = rSges(5,1) * t247 + rSges(5,2) * t301 + rSges(5,3) * t230;
t226 = rSges(6,2) * t230 + t160 * t241 + t247 * t288;
t225 = (-t196 / 0.2e1 + t195 / 0.2e1 + t194 / 0.2e1) * qJD(3);
t224 = t278 * t246;
t163 = t145 * t281 + t146 * t240 - t138;
t161 = t163 * t159;
t3 = t240 * t231 + (t161 - t253) * qJD(1) + t226 + t235;
t4 = (-t242 - t241 + (t146 * t281 + t308) * qJD(3)) * t159 + (t160 * t163 - t299) * qJD(1) + t234;
t223 = t159 * t3 + t160 * t4;
t222 = -t116 - t229;
t221 = -rSges(4,2) * t145 + t283;
t27 = (-t158 + t288) * t160 + t161;
t184 = t228 + t91;
t28 = t184 + t276;
t204 = t159 * t28 + t160 * t27;
t201 = t145 * t279 - t138;
t164 = (rSges(5,2) - pkin(3)) * t146 + t201;
t29 = (rSges(5,1) - t158) * t160 + t164 * t159;
t85 = rSges(5,3) * t256 + t298;
t30 = t184 + t85;
t203 = t159 * t30 + t160 * t29;
t51 = t222 * t159;
t52 = t222 * t160;
t202 = t159 * t51 + t160 * t52;
t83 = rSges(4,1) * t254 + t297;
t183 = -t138 - t221;
t176 = qJD(3) * t118;
t173 = qJD(3) * (Icges(5,4) * t145 + Icges(5,5) * t146);
t172 = qJD(3) * (Icges(6,4) * t146 - Icges(6,5) * t145);
t171 = qJD(3) * (-Icges(4,5) * t145 - Icges(4,6) * t146);
t5 = t159 * t275 + t160 * t276 + t284;
t165 = -qJ(5) * t245 - t218 * qJD(3) - qJD(5) * t145 - t81;
t92 = t116 * t248;
t8 = t160 * t165 + t229 * t248 + t92;
t9 = qJD(1) * t52 + t159 * t165;
t166 = qJD(3) * t5 + t159 * t9 + t160 * t8;
t162 = rSges(4,2) * t233 + rSges(4,3) * t247 - t160 * t176;
t60 = -t159 * t185 + t160 * t280;
t104 = t221 * qJD(3);
t87 = -rSges(5,1) * t160 + t159 * t220;
t82 = -rSges(4,3) * t160 + t159 * t221;
t58 = t252 * t159;
t57 = t228 + t83;
t56 = (rSges(4,3) - t158) * t160 + t183 * t159;
t55 = -qJD(1) * t61 + t148;
t54 = qJD(1) * t60 + t147;
t53 = t239 * t145 * t245;
t50 = qJD(1) * t294 + t160 * t173;
t49 = t159 * t173 + t257;
t48 = -qJD(1) * t77 + t160 * t172;
t47 = t159 * t172 + t258;
t34 = t159 * t171 + t259;
t33 = -qJD(1) * t62 + t160 * t171;
t26 = t118 * t244 + (t160 * t183 - t149) * qJD(1) + t250;
t25 = t147 + (-t253 + (-t138 - t283) * t159) * qJD(1) + t162;
t24 = qJD(1) * t59 + t159 * t277;
t23 = t160 * t277 - t219 * t248 + t92;
t22 = t159 * t205 + t160 * t294;
t21 = -t160 * t78 + t178;
t20 = t159 * t207 - t160 * t77;
t19 = -t160 * t76 + t180;
t18 = -t159 * t294 + t177;
t17 = t159 * t78 + t206 * t160;
t16 = t159 * t77 + t179;
t15 = t159 * t76 + t208 * t160;
t14 = t159 * t63 - t209 * t160;
t13 = t159 * t62 - t182;
t12 = -t160 * t63 - t181;
t11 = -t159 * t210 - t160 * t62;
t10 = t159 * t87 + t160 * t85 + t284;
t7 = (-t242 + (t146 * t279 - t282) * qJD(3)) * t159 + (t160 * t164 - t151) * qJD(1) + t234;
t6 = -pkin(3) * t231 + (-t253 + (t201 - t285) * t159) * qJD(1) + t227 + t235;
t2 = (qJD(1) * t87 + t227) * t160 + (t219 * t244 + (-t85 - t91 + t298) * qJD(1)) * t159 + t238;
t1 = (qJD(1) * t275 - t160 * t224 + t226) * t160 + ((rSges(6,2) * qJD(3) + qJD(5)) * t255 + (-t276 - t91 + t299) * qJD(1) - t159 * t224) * t159 + t238;
t31 = [0.2e1 * m(3) * (t54 * t61 + t55 * t60) + (t25 * t57 + t26 * t56) * t293 + (t29 * t7 + t30 * t6) * t292 + (t27 * t4 + t28 * t3) * t291 + (t200 + t193 + t186 - t309) * t246 + (Icges(6,3) * t145 + t189 - t190 + t198 - t261 + t310) * t245; m(3) * (t159 * t55 - t160 * t54 + (t159 * t61 + t160 * t60) * qJD(1)) + m(4) * (t159 * t26 - t160 * t25 + (t159 * t57 + t160 * t56) * qJD(1)) + m(5) * (qJD(1) * t203 + t159 * t7 - t160 * t6) + m(6) * (qJD(1) * t204 + t159 * t4 - t160 * t3); 0; m(5) * (t23 * t29 + t24 * t30 + t58 * t6 + t59 * t7) + m(6) * (t27 * t8 + t28 * t9 + t3 * t51 + t4 * t52) + (m(4) * (-t104 * t56 - t118 * t26) + t225 * t160 + (t65 * t303 + (t70 + t72) * t302) * t146 + t305 * t306) * t160 + (m(4) * (-t104 * t57 - t118 * t25) + t225 * t159 + (t302 * t73 + (t296 + t64) * t303) * t146 + t304 * t306) * t159 + ((t74 * t302 + t307 * t305 + (t68 + t67) * t303) * t160 + (t307 * t304 + (t69 + t295 + t66) * t303) * t159) * t145 + (-t179 / 0.2e1 - t177 / 0.2e1 + t182 / 0.2e1 + t180 / 0.2e1 + t178 / 0.2e1 - t181 / 0.2e1) * qJD(3) + ((-t57 * t287 + (t65 / 0.2e1 - t72 / 0.2e1 - t70 / 0.2e1) * t146 + (t67 / 0.2e1 + t68 / 0.2e1 - t74 / 0.2e1) * t145) * t160 + (t56 * t287 + (-t73 / 0.2e1 + t296 / 0.2e1 + t64 / 0.2e1) * t146 + (t69 / 0.2e1 + t295 / 0.2e1 + t66 / 0.2e1) * t145) * t159) * qJD(1); m(5) * (t23 * t159 - t160 * t24 + (t159 * t58 + t160 * t59) * qJD(1)) + m(6) * (qJD(1) * t202 + t8 * t159 - t160 * t9); (t10 * t2 + t23 * t59 + t24 * t58) * t292 + (t1 * t5 + t51 * t9 + t52 * t8) * t291 - t160 * ((t160 * t34 + (t12 + t182) * qJD(1)) * t160 + (t11 * qJD(1) + (-t245 * t65 - t246 * t67 + t259) * t159 + (-t33 + (t145 * t66 + t146 * t64) * qJD(3) - t209 * qJD(1)) * t160) * t159) + t159 * ((t159 * t50 + (t18 - t178) * qJD(1)) * t159 + (t17 * qJD(1) + (t245 * t296 + t246 * t295) * t160 + (-t49 + (t145 * t74 + t146 * t70) * qJD(3) + (t205 + t78) * qJD(1)) * t159) * t160) - t160 * ((t160 * t47 + (t19 - t179) * qJD(1)) * t160 + (t20 * qJD(1) + (t245 * t72 - t246 * t68 + t258) * t159 + (-t48 + (t145 * t69 - t146 * t73) * qJD(3) + t208 * qJD(1)) * t160) * t159) - t160 * ((t160 * t49 + (t21 - t177) * qJD(1)) * t160 + (t22 * qJD(1) + (t245 * t70 + t246 * t74 + t257) * t159 + (-t50 + (t145 * t295 + t146 * t296) * qJD(3) + t206 * qJD(1)) * t160) * t159) + t159 * ((t159 * t48 + (t16 - t180) * qJD(1)) * t159 + (t15 * qJD(1) + (-t245 * t73 + t246 * t69) * t160 + (-t47 + (-t145 * t68 + t146 * t72) * qJD(3) + (t207 + t76) * qJD(1)) * t159) * t160) + t159 * ((t159 * t33 + (t13 + t181) * qJD(1)) * t159 + (t14 * qJD(1) + (t245 * t64 + t246 * t66) * t160 + (-t34 + (-t145 * t67 - t146 * t65) * qJD(3) + (-t210 + t63) * qJD(1)) * t159) * t160) + ((t159 * t82 + t160 * t83) * ((qJD(1) * t82 + t162) * t160 + (-t159 * t176 + (-t83 + t297) * qJD(1)) * t159) + t249 * t118 * t104) * t293 + ((-t11 - t20 - t22) * t160 + (t12 + t19 + t21) * t159) * t248 + ((-t13 - t16 - t18) * t160 + (t14 + t15 + t17) * t159) * t247; 0.2e1 * (t203 * t290 + t204 * t289) * t245 + 0.2e1 * ((t159 * t6 + t160 * t7 + t247 * t30 - t248 * t29) * t290 + (t247 * t28 - t248 * t27 + t223) * t289) * t145; 0; 0.2e1 * ((t243 * t59 + t244 * t58 - t2) * t290 + (t243 * t52 + t244 * t51 - t1) * t289) * t146 + 0.2e1 * ((qJD(3) * t10 + t159 * t24 + t160 * t23 + t247 * t58 - t248 * t59) * t290 + (t247 * t51 - t248 * t52 + t166) * t289) * t145; 0.4e1 * (t290 + t289) * t53; m(6) * (-t204 * t246 + ((-t159 * t27 + t160 * t28) * qJD(1) + t223) * t146); 0; m(6) * ((-qJD(3) * t202 + t1) * t145 + ((-t159 * t52 + t160 * t51) * qJD(1) + t166) * t146); m(6) * (-t145 ^ 2 + t146 ^ 2) * t239 * qJD(3); -0.2e1 * m(6) * t53;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t31(1), t31(2), t31(4), t31(7), t31(11); t31(2), t31(3), t31(5), t31(8), t31(12); t31(4), t31(5), t31(6), t31(9), t31(13); t31(7), t31(8), t31(9), t31(10), t31(14); t31(11), t31(12), t31(13), t31(14), t31(15);];
Mq = res;

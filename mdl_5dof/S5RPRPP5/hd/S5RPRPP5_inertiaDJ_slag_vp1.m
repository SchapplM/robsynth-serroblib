% Calculate time derivative of joint inertia matrix for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:20
% DurationCPUTime: 7.95s
% Computational Cost: add. (1988->419), mult. (5132->575), div. (0->0), fcn. (3870->4), ass. (0->213)
t159 = sin(qJ(3));
t161 = cos(qJ(3));
t270 = Icges(5,5) * t161;
t273 = Icges(6,4) * t161;
t322 = t270 + t273 + (Icges(6,2) + Icges(5,3)) * t159;
t162 = cos(qJ(1));
t312 = rSges(6,3) + qJ(5);
t224 = t162 * t312;
t160 = sin(qJ(1));
t257 = t159 * t160;
t292 = rSges(6,1) + pkin(4);
t321 = -t292 * t257 + t224;
t271 = Icges(5,5) * t159;
t274 = Icges(6,4) * t159;
t277 = Icges(4,4) * t159;
t320 = t271 + t274 - t277 + (Icges(4,1) + Icges(5,1) + Icges(6,1)) * t161;
t319 = t322 * qJD(3);
t313 = qJD(1) / 0.2e1;
t318 = t161 * t313;
t317 = t320 * qJD(3);
t316 = t160 / 0.2e1;
t244 = qJD(3) * t162;
t315 = -t244 / 0.2e1;
t314 = -qJD(1) / 0.2e1;
t245 = qJD(3) * t161;
t248 = qJD(1) * t162;
t311 = t159 * t248 + t160 * t245;
t118 = pkin(3) * t161 + qJ(4) * t159;
t214 = rSges(4,1) * t159 + rSges(4,2) * t161;
t177 = t162 * t214;
t190 = Icges(4,2) * t161 + t277;
t71 = Icges(4,6) * t162 + t160 * t190;
t276 = Icges(4,4) * t161;
t196 = Icges(4,1) * t159 + t276;
t77 = Icges(4,5) * t162 + t160 * t196;
t201 = t159 * t77 + t161 * t71;
t172 = t201 * t162;
t186 = -Icges(6,2) * t161 + t274;
t67 = -Icges(6,6) * t162 + t160 * t186;
t192 = Icges(6,1) * t159 - t273;
t73 = -Icges(6,5) * t162 + t160 * t192;
t203 = t159 * t73 - t161 * t67;
t174 = t203 * t162;
t182 = -Icges(5,3) * t161 + t271;
t63 = Icges(5,6) * t162 + t160 * t182;
t194 = Icges(5,1) * t159 - t270;
t75 = Icges(5,4) * t162 + t160 * t194;
t205 = t159 * t75 - t161 * t63;
t176 = t205 * t162;
t254 = rSges(5,1) * t257 + rSges(5,2) * t162;
t284 = rSges(6,2) * t159;
t308 = t161 * t292 + t284;
t180 = Icges(6,5) * t159 - Icges(6,6) * t161;
t307 = Icges(6,3) * t160 + t162 * t180;
t306 = -Icges(5,6) * t160 + t162 * t182;
t184 = Icges(4,5) * t159 + Icges(4,6) * t161;
t305 = -Icges(4,3) * t160 + t162 * t184;
t304 = Icges(6,6) * t160 + t162 * t186;
t188 = Icges(5,4) * t159 - Icges(5,6) * t161;
t303 = -Icges(5,2) * t160 + t162 * t188;
t302 = -Icges(4,6) * t160 + t162 * t190;
t301 = Icges(6,5) * t160 + t162 * t192;
t300 = -Icges(5,4) * t160 + t162 * t194;
t299 = -Icges(4,5) * t160 + t162 * t196;
t298 = 2 * m(4);
t297 = 2 * m(5);
t296 = 2 * m(6);
t157 = t160 ^ 2;
t158 = t162 ^ 2;
t295 = m(5) / 0.2e1;
t294 = m(6) / 0.2e1;
t293 = -pkin(1) - pkin(6);
t291 = rSges(3,2) - pkin(1);
t285 = rSges(4,2) * t159;
t121 = rSges(4,1) * t161 - t285;
t290 = m(4) * t121;
t261 = qJ(4) * t161;
t88 = qJD(4) * t159 + (-pkin(3) * t159 + t261) * qJD(3);
t288 = t118 * t248 + t160 * t88;
t255 = t160 * t161;
t144 = pkin(3) * t257;
t239 = qJ(4) * t255;
t90 = t144 - t239;
t287 = rSges(5,3) * t255 - t254 - t90;
t286 = rSges(5,1) * t159;
t283 = rSges(6,2) * t161;
t282 = rSges(5,3) * t161;
t281 = t160 * rSges(5,2);
t280 = t160 * rSges(4,3);
t153 = t162 * rSges(4,3);
t212 = -rSges(6,1) * t159 + t283;
t256 = t159 * t162;
t278 = -pkin(4) * t256 - t160 * t312 + t162 * t212;
t61 = -Icges(6,3) * t162 + t160 * t180;
t260 = qJD(1) * t61;
t65 = Icges(4,3) * t162 + t160 * t184;
t259 = qJD(1) * t65;
t69 = Icges(5,2) * t162 + t160 * t188;
t258 = qJD(1) * t69;
t145 = pkin(3) * t256;
t152 = t162 * qJ(2);
t253 = t145 + t152;
t252 = qJ(2) * t248 + qJD(2) * t160;
t251 = pkin(1) * t162 + qJ(2) * t160;
t250 = t157 + t158;
t249 = qJD(1) * t160;
t247 = qJD(3) * t159;
t246 = qJD(3) * t160;
t243 = qJD(4) * t161;
t242 = -rSges(5,2) + t293;
t241 = -rSges(4,3) + t293;
t240 = rSges(6,2) * t255 + t321 - t90;
t231 = t159 * t246;
t238 = pkin(3) * t311 + qJ(4) * t231;
t237 = qJD(1) * t239 + t118 * t244;
t236 = -rSges(5,1) * t311 - rSges(5,3) * t231;
t229 = t161 * t248;
t235 = rSges(4,1) * t311 + rSges(4,2) * t229;
t81 = rSges(4,1) * t257 + rSges(4,2) * t255 + t153;
t234 = pkin(6) * t162 + t251;
t233 = t292 * t159;
t228 = t312 + t293;
t119 = rSges(6,1) * t161 + t284;
t227 = pkin(4) * t161 + t119;
t226 = (-rSges(6,2) - qJ(4)) * t161;
t225 = (-rSges(5,3) - qJ(4)) * t161;
t106 = t214 * qJD(3);
t223 = t106 * t250;
t222 = t242 * t160;
t150 = qJD(2) * t162;
t221 = t150 + t237;
t220 = t144 + t234;
t219 = (t180 / 0.2e1 - t188 / 0.2e1 - t184 / 0.2e1) * qJD(3);
t163 = -t160 * t243 + t238 + t252;
t164 = -rSges(6,2) * t231 - t249 * t312 - t292 * t311;
t3 = -qJD(5) * t162 + (t160 * t293 + t162 * t226) * qJD(1) + t163 - t164;
t4 = qJD(5) * t160 + (qJD(3) * t308 - t243) * t162 + (t228 * t162 + (t283 - qJ(2) + (-pkin(3) - t292) * t159) * t160) * qJD(1) + t221;
t218 = t160 * t4 - t162 * t3;
t6 = (t162 * t225 + t222) * qJD(1) + t163 - t236;
t120 = rSges(5,1) * t161 + rSges(5,3) * t159;
t7 = (qJD(3) * t120 - t243) * t162 + (t242 * t162 + (t282 - qJ(2) + (-rSges(5,1) - pkin(3)) * t159) * t160) * qJD(1) + t221;
t217 = t160 * t7 - t162 * t6;
t104 = t212 * qJD(3);
t10 = t119 * t248 + t160 * t104 + (t229 - t231) * pkin(4) + t288;
t215 = t227 * t160;
t92 = t118 * t249;
t9 = t92 + qJD(1) * t215 + (pkin(4) * t247 - t104 - t88) * t162;
t216 = t10 * t160 - t162 * t9;
t213 = t282 - t286;
t204 = -t159 * t300 + t161 * t306;
t202 = -t159 * t301 + t161 * t304;
t200 = -t159 * t299 - t161 * t302;
t105 = t213 * qJD(3);
t23 = t120 * t249 + t92 + (-t105 - t88) * t162;
t24 = t105 * t160 + t120 * t248 + t288;
t199 = t160 * t24 - t162 * t23;
t27 = (t226 + t233) * t162 + t228 * t160 + t253;
t28 = t160 * t226 + t220 - t321;
t198 = t160 * t27 - t162 * t28;
t191 = -Icges(4,2) * t159 + t276;
t189 = Icges(5,4) * t161 + Icges(5,6) * t159;
t185 = Icges(4,5) * t161 - Icges(4,6) * t159;
t181 = Icges(6,5) * t161 + Icges(6,6) * t159;
t179 = (t295 + t294) * t247;
t178 = rSges(3,3) * t162 + t160 * t291;
t175 = t204 * t160;
t173 = t202 * t160;
t171 = t200 * t160;
t167 = qJD(3) * t191;
t94 = t160 * t118;
t91 = t162 * t261 - t145;
t87 = t162 * t91;
t86 = -rSges(3,2) * t162 + rSges(3,3) * t160 + t251;
t85 = t152 + t178;
t84 = t280 - t177;
t83 = t162 * t213 + t281;
t59 = (-t118 - t120) * t162;
t58 = t120 * t160 + t94;
t56 = t150 + (t291 * t162 + (-rSges(3,3) - qJ(2)) * t160) * qJD(1);
t55 = qJD(1) * t178 + t252;
t54 = t234 + t81;
t53 = t160 * t241 + t152 + t177;
t52 = (-t118 - t227) * t162;
t51 = t94 + t215;
t42 = qJD(1) * t303 + t189 * t246;
t41 = -t189 * t244 + t258;
t38 = qJD(1) * t305 + t185 * t246;
t37 = -t185 * t244 + t259;
t34 = qJD(1) * t307 + t181 * t246;
t33 = -t181 * t244 + t260;
t32 = (-qJ(4) * t248 - qJD(4) * t160) * t161 + t238;
t31 = t162 * (qJD(1) * t144 + t162 * t243 - t237);
t30 = t160 * t225 + t220 + t254;
t29 = (t225 + t286) * t162 + t222 + t253;
t26 = t150 + t121 * t244 + (t241 * t162 + (-qJ(2) - t214) * t160) * qJD(1);
t25 = (-rSges(4,2) * t247 + qJD(1) * t241) * t160 + t235 + t252;
t22 = -t160 * t305 - t162 * t200;
t21 = t160 * t65 - t172;
t20 = -t160 * t303 - t162 * t204;
t19 = t160 * t69 - t176;
t18 = t160 * t307 - t162 * t202;
t17 = -t160 * t61 - t174;
t16 = -t162 * t305 + t171;
t15 = t160 * t201 + t162 * t65;
t14 = -t162 * t303 + t175;
t13 = t160 * t205 + t162 * t69;
t12 = t162 * t307 + t173;
t11 = t160 * t203 - t162 * t61;
t8 = t160 * t287 + t162 * t83 + t87;
t5 = t160 * t240 + t162 * t278 + t87;
t2 = t31 + (-t32 + (-t83 - t91 + t281) * qJD(1) + t236) * t160 + (-t120 * t244 + (t287 + t254) * qJD(1)) * t162;
t1 = t31 + (-t32 + (-t91 - t278) * qJD(1) + t164) * t160 + (-t308 * t244 + (t160 * t233 - t224 + t240) * qJD(1)) * t162;
t35 = [0.2e1 * m(3) * (t55 * t86 + t56 * t85) + (t25 * t54 + t26 * t53) * t298 + (t29 * t7 + t30 * t6) * t297 + (t27 * t4 + t28 * t3) * t296 + (t190 - t182 - t186) * t247 + (-t196 - t194 - t192) * t245 + (-t167 + t319) * t161 - t317 * t159; m(3) * (t160 * t56 - t162 * t55 + (t160 * t86 + t162 * t85) * qJD(1)) + m(4) * (t160 * t26 - t162 * t25 + (t160 * t54 + t162 * t53) * qJD(1)) + m(5) * ((t160 * t30 + t162 * t29) * qJD(1) + t217) + m(6) * ((t160 * t28 + t162 * t27) * qJD(1) + t218); 0; m(5) * (t23 * t30 + t24 * t29 + t58 * t7 + t59 * t6) + m(6) * (t10 * t27 + t28 * t9 + t3 * t52 + t4 * t51) + (m(4) * (t106 * t54 - t121 * t25) + t219 * t162 + t317 * t161 * t316 + (t299 + t300 + t301) * t318) * t162 + (m(4) * (-t106 * t53 + t121 * t26) + t219 * t160 + t320 * t161 * t315 + (t73 + t75 + t77) * t318) * t160 + ((t302 * t314 - t160 * t167 / 0.2e1 + t319 * t316) * t162 + (t71 * t314 + t191 * t244 / 0.2e1 + t322 * t315) * t160 + ((t304 + t306) * t162 + (t67 + t63) * t160) * t313) * t159 + (-t174 / 0.2e1 - t176 / 0.2e1 - t172 / 0.2e1 + (-t202 / 0.2e1 - t204 / 0.2e1 - t200 / 0.2e1) * t160) * qJD(3) + ((t54 * t290 + (-t73 / 0.2e1 - t75 / 0.2e1 - t77 / 0.2e1) * t161 + (-t67 / 0.2e1 - t63 / 0.2e1 + t71 / 0.2e1) * t159) * t160 + (t53 * t290 + (-t301 / 0.2e1 - t300 / 0.2e1 - t299 / 0.2e1) * t161 + (-t304 / 0.2e1 - t306 / 0.2e1 + t302 / 0.2e1) * t159) * t162) * qJD(1); m(5) * ((t160 * t59 + t162 * t58) * qJD(1) + t199) + m(6) * ((t160 * t52 + t162 * t51) * qJD(1) + t216) - m(4) * t223; (t1 * t5 + t10 * t51 + t52 * t9) * t296 + (t2 * t8 + t23 * t59 + t24 * t58) * t297 + ((-t160 * t81 + t162 * t84) * (-t160 * t235 + (-t121 * t158 + t157 * t285) * qJD(3) + ((-t81 + t153) * t162 + (t177 - t84 + t280) * t160) * qJD(1)) - t121 * t223) * t298 + t160 * ((-t160 * t33 + (-t17 + t173) * qJD(1)) * t160 + (t18 * qJD(1) + (-t245 * t73 - t247 * t67 - t260) * t162 + (-t34 + (t159 * t304 + t161 * t301) * qJD(3) + t203 * qJD(1)) * t160) * t162) + t160 * ((t160 * t41 + (-t19 + t175) * qJD(1)) * t160 + (t20 * qJD(1) + (-t245 * t75 - t247 * t63 + t258) * t162 + (t42 + (t159 * t306 + t161 * t300) * qJD(3) + t205 * qJD(1)) * t160) * t162) + t160 * ((t160 * t37 + (-t21 + t171) * qJD(1)) * t160 + (t22 * qJD(1) + (-t245 * t77 + t247 * t71 + t259) * t162 + (t38 + (-t159 * t302 + t161 * t299) * qJD(3) + t201 * qJD(1)) * t160) * t162) + t162 * ((-t162 * t34 + (t12 + t174) * qJD(1)) * t162 + (-t11 * qJD(1) + (-t245 * t301 - t247 * t304) * t160 + (-t33 + (t159 * t67 + t161 * t73) * qJD(3) + (t202 + t61) * qJD(1)) * t162) * t160) + t162 * ((t162 * t42 + (t14 + t176) * qJD(1)) * t162 + (-t13 * qJD(1) + (-t245 * t300 - t247 * t306) * t160 + (t41 + (t159 * t63 + t161 * t75) * qJD(3) + (t204 - t69) * qJD(1)) * t162) * t160) + t162 * ((t162 * t38 + (t16 + t172) * qJD(1)) * t162 + (-t15 * qJD(1) + (-t245 * t299 + t247 * t302) * t160 + (t37 + (-t159 * t71 + t161 * t77) * qJD(3) + (t200 - t65) * qJD(1)) * t162) * t160) + ((-t11 - t13 - t15) * t162 + (-t12 - t14 - t16) * t160) * t249 + ((t17 + t19 + t21) * t162 + (t18 + t20 + t22) * t160) * t248; 0.2e1 * ((t160 * t29 - t162 * t30) * t295 + t198 * t294) * t247 + 0.2e1 * ((-t248 * t29 - t249 * t30 - t217) * t295 + (-t248 * t27 - t249 * t28 - t218) * t294) * t161; 0.2e1 * t250 * t179; 0.2e1 * ((-t244 * t52 + t246 * t51 + t1) * t294 + (-t244 * t59 + t246 * t58 + t2) * t295) * t159 + 0.2e1 * ((qJD(3) * t5 - t248 * t51 - t249 * t52 - t216) * t294 + (qJD(3) * t8 - t248 * t58 - t249 * t59 - t199) * t295) * t161; 0.4e1 * (0.1e1 - t250) * t161 * t179; m(6) * (qJD(1) * t198 - t160 * t3 - t162 * t4); 0; m(6) * (-t10 * t162 - t160 * t9 + (t160 * t51 - t162 * t52) * qJD(1)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t35(1), t35(2), t35(4), t35(7), t35(11); t35(2), t35(3), t35(5), t35(8), t35(12); t35(4), t35(5), t35(6), t35(9), t35(13); t35(7), t35(8), t35(9), t35(10), t35(14); t35(11), t35(12), t35(13), t35(14), t35(15);];
Mq = res;

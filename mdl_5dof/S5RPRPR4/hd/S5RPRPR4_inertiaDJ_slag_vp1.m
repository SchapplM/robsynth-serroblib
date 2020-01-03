% Calculate time derivative of joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:16
% EndTime: 2020-01-03 11:38:40
% DurationCPUTime: 7.12s
% Computational Cost: add. (7630->407), mult. (5952->544), div. (0->0), fcn. (4556->10), ass. (0->238)
t355 = Icges(4,3) + Icges(5,3);
t172 = qJ(3) + pkin(9);
t162 = sin(t172);
t164 = cos(t172);
t175 = sin(qJ(3));
t177 = cos(qJ(3));
t354 = Icges(4,5) * t177 + Icges(5,5) * t164 - Icges(4,6) * t175 - Icges(5,6) * t162;
t173 = qJ(1) + pkin(8);
t163 = sin(t173);
t165 = cos(t173);
t289 = Icges(4,4) * t175;
t216 = Icges(4,1) * t177 - t289;
t318 = Icges(4,5) * t163 + t165 * t216;
t287 = Icges(5,4) * t162;
t214 = Icges(5,1) * t164 - t287;
t319 = Icges(5,5) * t163 + t165 * t214;
t288 = Icges(4,4) * t177;
t211 = -Icges(4,2) * t175 + t288;
t321 = Icges(4,6) * t163 + t165 * t211;
t286 = Icges(5,4) * t164;
t209 = -Icges(5,2) * t162 + t286;
t322 = Icges(5,6) * t163 + t165 * t209;
t350 = -t162 * t322 + t164 * t319 - t175 * t321 + t177 * t318;
t353 = t163 * t350;
t101 = -Icges(4,6) * t165 + t163 * t211;
t103 = -Icges(4,5) * t165 + t163 * t216;
t91 = -Icges(5,6) * t165 + t163 * t209;
t93 = -Icges(5,5) * t165 + t163 * t214;
t349 = t101 * t175 - t103 * t177 + t162 * t91 - t164 * t93;
t352 = t349 * t163;
t351 = t354 * t163 - t355 * t165;
t347 = t355 * t163 + t354 * t165;
t348 = Icges(4,5) * t175 + Icges(5,5) * t162 + Icges(4,6) * t177 + Icges(5,6) * t164;
t166 = qJ(5) + t172;
t156 = sin(t166);
t301 = rSges(6,2) * t156;
t157 = cos(t166);
t304 = rSges(6,1) * t157;
t226 = -t301 + t304;
t87 = -rSges(6,3) * t165 + t163 * t226;
t202 = Icges(6,5) * t157 - Icges(6,6) * t156;
t326 = Icges(6,3) * t163 + t165 * t202;
t346 = t326 * qJD(1);
t345 = t165 * t347 - t353;
t343 = t349 * t165;
t342 = t351 * qJD(1);
t169 = cos(qJ(1)) * pkin(1);
t159 = qJD(1) * t169;
t160 = t163 ^ 2;
t315 = 2 * m(5);
t127 = rSges(5,1) * t162 + rSges(5,2) * t164;
t190 = qJD(3) * t127;
t155 = t165 * pkin(2);
t174 = -qJ(4) - pkin(6);
t168 = t177 * pkin(3);
t158 = t168 + pkin(2);
t253 = qJD(1) * t165;
t130 = t158 * t253;
t248 = qJD(3) * t175;
t244 = pkin(3) * t248;
t198 = -t163 * t244 + t130;
t254 = qJD(1) * t163;
t148 = t165 * t174;
t259 = t163 * t158 + t148;
t79 = -pkin(2) * t163 + pkin(6) * t165 + t259;
t138 = t165 * t158;
t154 = t163 * pkin(6);
t258 = t154 + t155;
t80 = t163 * t174 - t138 + t258;
t245 = t163 * (-qJD(4) * t165 + (-t155 + (-pkin(6) - t174) * t163) * qJD(1) + t198) + t80 * t254 + t79 * t253;
t306 = rSges(5,1) * t164;
t260 = rSges(5,3) * t254 + t253 * t306;
t328 = rSges(5,3) * t165 - t163 * t306;
t152 = pkin(6) * t253;
t153 = qJD(4) * t163;
t58 = t165 * t244 + t152 - t153 + (t148 + (-pkin(2) + t158) * t163) * qJD(1);
t302 = rSges(5,2) * t162;
t95 = -t163 * t302 - t328;
t227 = -t302 + t306;
t332 = t165 * t227;
t96 = -rSges(5,3) * t163 - t332;
t4 = (qJD(1) * t96 - t163 * t190 + t260) * t163 + (-t58 - t165 * t190 + (t95 + t328) * qJD(1)) * t165 + t245;
t341 = t315 * t4;
t340 = t351 * t165 + t352;
t339 = t347 * t163 + t350 * t165;
t338 = -t351 * t163 + t343;
t337 = -t348 * qJD(3) * t165 - t342;
t251 = qJD(3) * t163;
t336 = -t347 * qJD(1) + t348 * t251;
t303 = rSges(4,2) * t175;
t307 = rSges(4,1) * t177;
t145 = t163 * t307;
t299 = rSges(4,3) * t165;
t327 = t299 - t145;
t105 = -t163 * t303 - t327;
t228 = -t303 + t307;
t106 = -rSges(4,3) * t163 - t228 * t165;
t144 = rSges(4,1) * t175 + rSges(4,2) * t177;
t191 = qJD(3) * t144;
t179 = rSges(4,3) * t254 - t163 * t191 + t253 * t307;
t182 = t165 * t191;
t13 = (qJD(1) * t106 + t179) * t163 + (-t182 + (t105 + t327) * qJD(1)) * t165;
t316 = 2 * m(4);
t335 = t13 * t316;
t285 = Icges(6,4) * t156;
t212 = Icges(6,1) * t157 - t285;
t320 = Icges(6,5) * t163 + t165 * t212;
t284 = Icges(6,4) * t157;
t207 = -Icges(6,2) * t156 + t284;
t323 = Icges(6,6) * t163 + t165 * t207;
t224 = -t156 * t323 + t157 * t320;
t334 = t163 * t224;
t170 = -pkin(7) + t174;
t256 = t170 - t174;
t333 = t163 * t256;
t171 = qJD(3) + qJD(5);
t108 = t207 * t171;
t109 = t212 * t171;
t118 = Icges(6,5) * t156 + Icges(6,6) * t157;
t119 = Icges(6,2) * t157 + t285;
t120 = Icges(6,1) * t156 + t284;
t317 = -qJD(1) * t118 + (t119 * t171 - t109) * t157 + (t120 * t171 + t108) * t156;
t314 = 2 * m(6);
t161 = t165 ^ 2;
t313 = -t163 / 0.2e1;
t312 = -t165 / 0.2e1;
t311 = m(4) * t144;
t310 = pkin(3) * t175;
t309 = pkin(4) * t162;
t308 = pkin(4) * t164;
t167 = sin(qJ(1)) * pkin(1);
t305 = rSges(6,1) * t156;
t296 = t162 * t93;
t295 = t162 * t319;
t121 = rSges(6,2) * t157 + t305;
t231 = -t309 - t310;
t197 = -t121 + t231;
t69 = t197 * t163;
t294 = t163 * t69;
t293 = t164 * t91;
t292 = t164 * t322;
t291 = rSges(5,3) - t174;
t290 = rSges(6,3) - t170;
t81 = -Icges(6,3) * t165 + t163 * t202;
t274 = qJD(1) * t81;
t271 = t101 * t177;
t270 = t177 * t321;
t269 = t103 * t175;
t268 = t175 * t318;
t267 = t156 * t171;
t266 = t157 * t171;
t265 = t163 * t171;
t264 = t165 * t171;
t125 = t231 * qJD(3);
t134 = t158 + t308;
t263 = t163 * t125 + t134 * t253;
t262 = t163 * t134 + t165 * t170;
t261 = rSges(6,3) * t254 + t253 * t304;
t257 = t160 + t161;
t252 = qJD(3) * t162;
t250 = qJD(3) * t164;
t247 = qJD(3) * t177;
t88 = -rSges(6,3) * t163 - t165 * t226;
t246 = t163 * (-t265 * t305 + (-t156 * t253 - t157 * t265) * rSges(6,2) + t261) + t88 * t254 + t87 * t253;
t83 = -Icges(6,6) * t165 + t163 * t207;
t85 = -Icges(6,5) * t165 + t163 * t212;
t225 = t156 * t83 - t157 * t85;
t16 = -t163 * t225 - t165 * t81;
t17 = t165 * t326 - t334;
t193 = t225 * t165;
t18 = -t163 * t81 + t193;
t19 = t163 * t326 + t165 * t224;
t44 = qJD(1) * t85 + t120 * t264;
t235 = t171 * t323 + t44;
t42 = qJD(1) * t83 + t119 * t264;
t237 = -t171 * t320 + t42;
t40 = t118 * t264 + t274;
t41 = -t118 * t265 + t346;
t43 = t323 * qJD(1) - t119 * t265;
t45 = t320 * qJD(1) - t120 * t265;
t243 = -t163 * ((t163 * t40 + (t18 + t334) * qJD(1)) * t163 + (-t19 * qJD(1) + (-t156 * t43 + t157 * t45 - t266 * t83 - t267 * t85 + t274) * t165 + (t41 + t235 * t157 - t237 * t156 + (t225 - t326) * qJD(1)) * t163) * t165) + (-t16 * t165 - t163 * t17) * t254;
t242 = pkin(2) - t303;
t234 = t171 * t83 - t45;
t236 = t171 * t85 + t43;
t2 = (t165 * t41 + (-t17 + t193) * qJD(1)) * t165 + (t16 * qJD(1) + (t156 * t42 - t157 * t44 - t266 * t323 - t267 * t320 + t346) * t163 + (t40 + t234 * t157 + t236 * t156 + (t224 - t81) * qJD(1)) * t165) * t163;
t6 = -t163 * t19 - t165 * t18;
t241 = -qJD(1) * t6 - t2;
t238 = -t127 - t310;
t97 = t238 * t163;
t49 = t167 + t87 + t262;
t196 = t134 + t226;
t50 = t163 * t290 + t165 * t196 + t169;
t219 = t163 * t50 - t165 * t49;
t215 = Icges(4,1) * t175 + t288;
t213 = Icges(5,1) * t162 + t286;
t210 = Icges(4,2) * t177 + t289;
t208 = Icges(5,2) * t164 + t287;
t199 = t119 * t156 - t120 * t157;
t194 = t121 * t171;
t181 = -qJD(1) * t199 - t202 * t171;
t189 = (t156 * t235 + t157 * t237 + t181 * t163 + t317 * t165) * t313 + (-t156 * t234 + t157 * t236 - t317 * t163 + t181 * t165) * t312 + (-t118 * t165 + t156 * t85 + t157 * t83 - t163 * t199) * t254 / 0.2e1 - (-t118 * t163 - t156 * t320 - t157 * t323 + t165 * t199) * t253 / 0.2e1;
t187 = qJD(3) * t215;
t186 = qJD(3) * t213;
t185 = qJD(3) * t210;
t184 = qJD(3) * t208;
t183 = qJD(1) * t197;
t180 = t238 * qJD(3);
t147 = t165 * t310;
t140 = t165 * pkin(3) * t247;
t135 = t228 * qJD(3);
t117 = t227 * qJD(3);
t110 = t226 * t171;
t98 = t127 * t165 + t147;
t78 = t163 * t87;
t77 = t163 * t79;
t72 = -t106 + t169 + t258;
t71 = t145 + t167 + (-rSges(4,3) - pkin(6)) * t165 + t242 * t163;
t70 = t147 + (t121 + t309) * t165;
t68 = -t134 * t165 + t138 + t333;
t67 = -t259 + t262;
t60 = t163 * t291 + t138 + t169 + t332;
t59 = t167 + t95 + t259;
t48 = -t127 * t253 - t117 * t163 + (-t163 * t247 - t175 * t253) * pkin(3);
t47 = qJD(1) * t97 + t117 * t165 + t140;
t46 = t87 * qJD(1) + t121 * t264;
t39 = t159 + (t165 * t242 + t154) * qJD(1) + t179;
t38 = t152 - t182 + (t299 - t167 + (-pkin(2) - t228) * t163) * qJD(1);
t36 = -t165 * t88 + t78;
t31 = (-t110 + (-t168 - t308) * qJD(3)) * t163 + t165 * t183;
t30 = t140 + (pkin(4) * t250 + t110) * t165 + t163 * t183;
t25 = t130 + t159 + (-qJD(1) * t302 - qJD(4)) * t165 + (-qJD(1) * t174 + t180) * t163 + t260;
t24 = t153 + t165 * t180 + (-t167 + t291 * t165 + (-t158 - t227) * t163) * qJD(1);
t15 = t159 + (-qJD(1) * t301 - qJD(4)) * t165 + (-qJD(1) * t170 - t194) * t163 + t261 + t263;
t14 = t153 + (t125 - t194) * t165 + (-t163 * t196 + t165 * t290 - t167) * qJD(1);
t12 = -t165 * t46 + t246;
t7 = t163 * t67 + t77 + t78 + (-t68 - t80 - t88) * t165;
t3 = ((t68 - t333) * qJD(1) - t198 + t263) * t163 + (-t46 - t58 + (t125 + t244) * t165 + (t67 - t256 * t165 + (-t134 + t158) * t163) * qJD(1)) * t165 + t245 + t246;
t1 = [t177 * t187 + t216 * t248 - t175 * t185 + t211 * t247 + t164 * t186 + t214 * t252 - t162 * t184 + t209 * t250 + (t38 * t72 + t39 * t71) * t316 + (t24 * t60 + t59 * t25) * t315 + t120 * t266 + t156 * t109 - t119 * t267 + t157 * t108 + (t14 * t50 + t15 * t49) * t314; 0; 0; ((t296 / 0.2e1 + t293 / 0.2e1 - t71 * t311 + t271 / 0.2e1 + t269 / 0.2e1) * t163 + (-t72 * t311 + t295 / 0.2e1 + t292 / 0.2e1 + t270 / 0.2e1 + t268 / 0.2e1) * t165) * qJD(1) + m(5) * (t24 * t97 + t25 * t98 + t47 * t59 + t48 * t60) + m(6) * (t14 * t69 + t15 * t70 + t30 * t49 + t31 * t50) + m(4) * ((-t163 * t38 + t165 * t39) * t144 + (-t163 * t72 + t165 * t71) * t135) + t189 + (-t350 * qJD(3) + t162 * (qJD(1) * t93 + t165 * t186) + t164 * (qJD(1) * t91 + t165 * t184) + t175 * (qJD(1) * t103 + t165 * t187) + t177 * (qJD(1) * t101 + t165 * t185)) * t313 + (-t349 * qJD(3) + t162 * (t319 * qJD(1) - t213 * t251) + t164 * (t322 * qJD(1) - t208 * t251) + t175 * (t318 * qJD(1) - t215 * t251) + t177 * (t321 * qJD(1) - t210 * t251)) * t312 + t354 * qJD(3) * (t161 / 0.2e1 + t160 / 0.2e1); m(4) * t13 + m(5) * t4 + m(6) * t3; (t7 * t3 + t30 * t70 + t31 * t69) * t314 - t6 * t253 + (t4 * t77 + t98 * t47 + t97 * t48) * t315 + t257 * t144 * t135 * t316 + t243 + (t105 * t335 + t337 * t160 + t339 * t253 + t95 * t341 + (-t338 - t345 - t353) * t254) * t163 + (-t2 + (-t80 - t96) * t341 - t106 * t335 + (t336 * t165 + (-t343 + t345) * qJD(1)) * t165 + t340 * t254 + t338 * t253 + ((t101 * t247 + t103 * t248 + t91 * t250 + t93 * t252 + t337 - t342) * t165 + (t247 * t321 + t248 * t318 + t250 * t322 + t252 * t319 + t336) * t163 + ((-t293 - t296 - t269 - t271) * t165 + (-t268 - t270 - t292 - t295) * t163) * qJD(3) + ((-t350 + t351) * t165 - t352 + t339 + t340) * qJD(1)) * t163) * t165; m(5) * (-t163 * t25 - t165 * t24 + (t163 * t60 - t165 * t59) * qJD(1)) + m(6) * (qJD(1) * t219 - t14 * t165 - t15 * t163); 0; m(6) * (-t163 * t30 - t165 * t31 + (-t165 * t70 + t294) * qJD(1)) + m(5) * (-t163 * t47 - t165 * t48 + (t163 * t97 - t165 * t98) * qJD(1)); 0; m(6) * (-t219 * t110 + (-t14 * t163 + t15 * t165 + (-t163 * t49 - t165 * t50) * qJD(1)) * t121) + t189; m(6) * t12; m(6) * (-t110 * t294 + t12 * t7 + t36 * t3 + (-t163 * t31 - t254 * t70) * t121) + (m(6) * (t110 * t70 + (-qJD(1) * t69 + t30) * t121) + t241) * t165 + t243; 0; (t110 * t121 * t257 + t36 * t12) * t314 + t241 * t165 + t243;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

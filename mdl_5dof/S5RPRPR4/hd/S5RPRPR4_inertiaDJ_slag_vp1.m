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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:53:10
% EndTime: 2019-12-05 17:53:25
% DurationCPUTime: 5.70s
% Computational Cost: add. (7630->404), mult. (5952->536), div. (0->0), fcn. (4556->10), ass. (0->241)
t355 = Icges(4,3) + Icges(5,3);
t171 = qJ(3) + pkin(9);
t163 = sin(t171);
t165 = cos(t171);
t174 = sin(qJ(3));
t176 = cos(qJ(3));
t354 = Icges(4,5) * t176 + Icges(5,5) * t165 - Icges(4,6) * t174 - Icges(5,6) * t163;
t172 = qJ(1) + pkin(8);
t164 = sin(t172);
t166 = cos(t172);
t353 = -t164 * t354 + t166 * t355;
t290 = Icges(5,4) * t165;
t201 = -Icges(5,2) * t163 + t290;
t87 = Icges(5,6) * t166 - t164 * t201;
t291 = Icges(5,4) * t163;
t206 = Icges(5,1) * t165 - t291;
t89 = Icges(5,5) * t166 - t164 * t206;
t292 = Icges(4,4) * t176;
t203 = -Icges(4,2) * t174 + t292;
t97 = Icges(4,6) * t166 - t164 * t203;
t293 = Icges(4,4) * t174;
t208 = Icges(4,1) * t176 - t293;
t99 = Icges(4,5) * t166 - t164 * t208;
t330 = t163 * t87 - t165 * t89 + t174 * t97 - t176 * t99;
t352 = t166 * t330;
t100 = Icges(4,5) * t164 + t166 * t208;
t88 = Icges(5,6) * t164 + t166 * t201;
t90 = Icges(5,5) * t164 + t166 * t206;
t98 = Icges(4,6) * t164 + t166 * t203;
t331 = t100 * t176 - t163 * t88 + t165 * t90 - t174 * t98;
t350 = t164 * t355 + t354 * t166;
t349 = t330 * t164;
t348 = Icges(4,5) * t174 + Icges(5,5) * t163 + Icges(4,6) * t176 + Icges(5,6) * t165;
t347 = -t164 * t353 + t352;
t346 = t331 * t164;
t252 = qJD(3) * t176;
t259 = qJD(1) * t166;
t345 = t164 * t252 + t174 * t259;
t255 = qJD(3) * t165;
t344 = t163 * t259 + t164 * t255;
t343 = t353 * qJD(1);
t167 = qJ(5) + t171;
t157 = sin(t167);
t158 = cos(t167);
t114 = rSges(6,1) * t157 + rSges(6,2) * t158;
t162 = t166 ^ 2;
t326 = 2 * m(5);
t256 = qJD(3) * t164;
t311 = rSges(5,1) * t163;
t248 = -rSges(5,2) * t344 - t256 * t311;
t155 = t166 * rSges(5,3);
t307 = rSges(5,2) * t163;
t266 = t164 * t307 + t155;
t301 = t164 * rSges(5,3);
t310 = rSges(5,1) * t165;
t173 = -qJ(4) - pkin(6);
t168 = t176 * pkin(3);
t159 = t168 + pkin(2);
t314 = pkin(2) - t159;
t75 = (-pkin(6) - t173) * t166 + t314 * t164;
t313 = t164 * t310 - t266 - t75;
t122 = rSges(5,2) * t165 + t311;
t254 = qJD(3) * t166;
t334 = t122 * t254;
t260 = qJD(1) * t164;
t152 = pkin(2) * t260;
t253 = qJD(3) * t174;
t258 = qJD(1) * t173;
t247 = t159 * t260 + (pkin(3) * t253 + t258) * t166;
t251 = qJD(4) * t164;
t191 = t247 - t251;
t51 = t166 * (-pkin(6) * t259 + t152 - t191);
t153 = qJD(4) * t166;
t241 = t164 * t253;
t265 = pkin(3) * t241 + t164 * t258;
t245 = t153 + t265;
t329 = pkin(6) * t164 + t166 * t314;
t58 = qJD(1) * t329 + t245;
t151 = t164 * t173;
t76 = -t151 - t329;
t222 = -t307 + t310;
t92 = t166 * t222 + t301;
t4 = t51 + (-t58 + (-t76 - t92 + t301) * qJD(1) + t248) * t164 + (-t334 + (t313 + t266) * qJD(1)) * t166;
t342 = t326 * t4;
t341 = -t166 * t353 - t349;
t340 = -t166 * t350 + t346;
t339 = t164 * t350 + t166 * t331;
t338 = -t254 * t348 + t343;
t337 = -qJD(1) * t350 + t256 * t348;
t271 = t164 * t174;
t263 = rSges(4,2) * t271 + rSges(4,3) * t166;
t312 = rSges(4,1) * t176;
t101 = -t164 * t312 + t263;
t223 = -rSges(4,2) * t174 + t312;
t302 = t164 * rSges(4,3);
t102 = t166 * t223 + t302;
t246 = -rSges(4,1) * t241 - rSges(4,2) * t345;
t148 = rSges(4,1) * t174 + rSges(4,2) * t176;
t333 = t148 * t254;
t13 = ((-t102 + t302) * qJD(1) + t246) * t164 + (-t333 + (-t101 + t263) * qJD(1)) * t166;
t327 = 2 * m(4);
t336 = t13 * t327;
t194 = Icges(6,5) * t158 - Icges(6,6) * t157;
t78 = Icges(6,3) * t164 + t166 * t194;
t335 = qJD(1) * t78;
t288 = Icges(6,4) * t158;
t199 = -Icges(6,2) * t157 + t288;
t79 = Icges(6,6) * t166 - t164 * t199;
t289 = Icges(6,4) * t157;
t204 = Icges(6,1) * t158 - t289;
t81 = Icges(6,5) * t166 - t164 * t204;
t220 = t157 * t79 - t158 * t81;
t332 = t166 * t220;
t320 = sin(qJ(1)) * pkin(1);
t160 = qJD(1) * t320;
t80 = Icges(6,6) * t164 + t166 * t199;
t82 = Icges(6,5) * t164 + t166 * t204;
t170 = qJD(3) + qJD(5);
t104 = t199 * t170;
t105 = t204 * t170;
t111 = Icges(6,5) * t157 + Icges(6,6) * t158;
t112 = Icges(6,2) * t158 + t289;
t113 = Icges(6,1) * t157 + t288;
t328 = t157 * (t113 * t170 + t104) + (t112 * t170 - t105) * t158 - qJD(1) * t111;
t325 = 2 * m(6);
t161 = t164 ^ 2;
t324 = t164 / 0.2e1;
t323 = t166 / 0.2e1;
t322 = -rSges(4,3) - pkin(6);
t321 = m(4) * t148;
t319 = cos(qJ(1)) * pkin(1);
t318 = pkin(3) * t174;
t317 = pkin(4) * t163;
t316 = pkin(4) * t165;
t308 = rSges(6,1) * t158;
t306 = rSges(6,2) * t157;
t304 = t163 * t89;
t303 = t163 * t90;
t300 = t164 * rSges(6,3);
t299 = t165 * t87;
t298 = t165 * t88;
t154 = t166 * rSges(6,3);
t297 = t174 * t99;
t296 = t176 * t97;
t295 = t176 * t98;
t169 = -pkin(7) + t173;
t294 = -rSges(6,3) + t169;
t77 = Icges(6,3) * t166 - t164 * t194;
t278 = qJD(1) * t77;
t275 = t100 * t174;
t274 = t157 * t170;
t273 = t158 * t170;
t272 = t164 * t170;
t270 = t166 * t169;
t269 = t166 * t170;
t268 = t164 * t306 + t154;
t131 = t159 + t316;
t267 = t131 - t159;
t264 = t345 * pkin(3);
t262 = t161 + t162;
t257 = qJD(3) * t163;
t83 = -t164 * t308 + t268;
t250 = -(-t169 + t173) * t166 + t267 * t164 - t75 - t83;
t150 = pkin(3) * t271;
t249 = t114 * t272 + t259 * t306;
t239 = -pkin(2) - t312;
t42 = qJD(1) * t79 - t112 * t269;
t236 = t170 * t82 + t42;
t43 = -qJD(1) * t80 + t112 * t272;
t235 = t170 * t81 + t43;
t44 = qJD(1) * t81 - t113 * t269;
t234 = -t170 * t80 + t44;
t45 = -qJD(1) * t82 + t113 * t272;
t233 = t170 * t79 - t45;
t232 = -t159 - t310;
t231 = -t131 - t308;
t227 = -t317 - t318;
t120 = t227 * qJD(3);
t228 = -t164 * t120 + t169 * t260;
t226 = -t301 - t319;
t225 = (t114 + t317) * t164;
t221 = -t306 + t308;
t219 = t157 * t80 - t158 * t82;
t49 = t164 * t231 + t268 - t270 - t320;
t190 = t131 + t221;
t50 = t164 * t294 - t166 * t190 - t319;
t214 = t164 * t50 - t166 * t49;
t69 = t150 + t225;
t70 = (-t114 + t227) * t166;
t213 = t164 * t69 - t166 * t70;
t207 = Icges(4,1) * t174 + t292;
t205 = Icges(5,1) * t163 + t290;
t202 = Icges(4,2) * t176 + t293;
t200 = Icges(5,2) * t165 + t291;
t193 = t112 * t157 - t113 * t158;
t16 = t164 * t220 + t166 * t77;
t187 = t219 * t164;
t17 = t166 * t78 + t187;
t18 = t164 * t77 - t332;
t19 = t164 * t78 - t166 * t219;
t40 = -t111 * t269 + t278;
t41 = t111 * t272 - t335;
t192 = -t260 * (t16 * t166 + t164 * t17) + t164 * ((t164 * t40 + (-t18 + t187) * qJD(1)) * t164 + (t19 * qJD(1) + (-t157 * t43 + t158 * t45 - t273 * t79 - t274 * t81 + t278) * t166 + (t41 + t234 * t158 - t236 * t157 + (t220 + t78) * qJD(1)) * t164) * t166) + t166 * ((t166 * t41 + (t17 + t332) * qJD(1)) * t166 + (-t16 * qJD(1) + (t157 * t42 - t158 * t44 + t273 * t80 + t274 * t82 - t335) * t164 + (t40 + t233 * t158 + t235 * t157 + (t219 - t77) * qJD(1)) * t166) * t164) + (t164 * t19 + t166 * t18) * t259;
t188 = t164 * t322 - t319;
t178 = qJD(1) * t193 + t170 * t194;
t184 = (t157 * t234 + t158 * t236 + t164 * t178 - t166 * t328) * t324 + (-t157 * t233 + t158 * t235 + t164 * t328 + t166 * t178) * t323 - (t111 * t166 + t157 * t81 + t158 * t79 + t164 * t193) * t260 / 0.2e1 + (t111 * t164 + t157 * t82 + t158 * t80 - t166 * t193) * t259 / 0.2e1;
t183 = qJD(3) * t207;
t182 = qJD(3) * t205;
t181 = qJD(3) * t202;
t180 = qJD(3) * t200;
t141 = qJD(1) * t150;
t132 = t223 * qJD(3);
t110 = t222 * qJD(3);
t106 = t221 * t170;
t94 = (-t122 - t318) * t166;
t93 = t122 * t164 + t150;
t84 = t166 * t221 + t300;
t74 = t166 * t84;
t73 = t166 * t76;
t72 = (-pkin(2) - t223) * t166 + t188;
t71 = pkin(6) * t166 + t164 * t239 + t263 - t320;
t68 = -t164 * t169 + t166 * t267 + t151;
t60 = t151 + (-t159 - t222) * t166 + t226;
t59 = t164 * t232 - t166 * t173 + t266 - t320;
t48 = t110 * t164 + t122 * t259 + t264;
t47 = t122 * t260 + t141 + (-pkin(3) * t252 - t110) * t166;
t46 = (-t166 * t308 - t300) * qJD(1) + t249;
t39 = (t166 * t239 + t188) * qJD(1) - t246;
t38 = t152 + t160 + t333 + (t164 * t223 + t166 * t322) * qJD(1);
t37 = t166 * (-t114 * t269 + (-t164 * t221 + t154) * qJD(1));
t36 = -t164 * t83 + t74;
t31 = pkin(4) * t344 + t106 * t164 + t114 * t259 + t264;
t30 = t141 + qJD(1) * t225 + (-t106 + (-t168 - t316) * qJD(3)) * t166;
t25 = (t166 * t232 + t226) * qJD(1) + t245 - t248;
t24 = t160 + t334 + (t164 * t222 - t155) * qJD(1) + t191;
t15 = t153 + (t166 * t231 - t300 - t319) * qJD(1) + t228 + t249;
t14 = -t251 + t160 + (t114 * t170 - t120) * t166 + (t164 * t190 + t166 * t294) * qJD(1);
t12 = -t164 * t46 + t37 + (-t164 * t84 - t166 * t83) * qJD(1);
t7 = t164 * t250 + t166 * t68 + t73 + t74;
t3 = t51 + t37 + (t166 * t120 + (t250 - t270) * qJD(1) + t247) * t166 + (-t46 - t58 + (-t159 * t166 - t68 - t76 - t84) * qJD(1) - t228 + t265) * t164;
t1 = [t176 * t183 + t208 * t253 - t174 * t181 + t203 * t252 + t165 * t182 + t206 * t257 - t163 * t180 + t201 * t255 + (t38 * t72 + t39 * t71) * t327 + (t24 * t60 + t59 * t25) * t326 + t113 * t273 + t157 * t105 - t112 * t274 + t158 * t104 + (t14 * t50 + t15 * t49) * t325; 0; 0; m(4) * ((t164 * t38 - t166 * t39) * t148 + (t164 * t72 - t166 * t71) * t132) + m(6) * (t14 * t69 + t15 * t70 + t30 * t49 + t31 * t50) + m(5) * (t24 * t93 + t25 * t94 + t47 * t59 + t48 * t60) + ((-t304 / 0.2e1 - t299 / 0.2e1 + t71 * t321 - t297 / 0.2e1 - t296 / 0.2e1) * t164 + (t72 * t321 + t303 / 0.2e1 + t298 / 0.2e1 + t275 / 0.2e1 + t295 / 0.2e1) * t166) * qJD(1) + t184 + (qJD(3) * t331 + t163 * (qJD(1) * t89 - t205 * t254) + t165 * (qJD(1) * t87 - t200 * t254) + t174 * (qJD(1) * t99 - t207 * t254) + t176 * (qJD(1) * t97 - t202 * t254)) * t324 + (-qJD(3) * t330 + t163 * (-qJD(1) * t90 + t164 * t182) + t165 * (-qJD(1) * t88 + t164 * t180) + t174 * (-qJD(1) * t100 + t164 * t183) + t176 * (-qJD(1) * t98 + t164 * t181)) * t323 + t354 * qJD(3) * (t162 / 0.2e1 + t161 / 0.2e1); m(4) * t13 + m(5) * t4 + m(6) * t3; (t7 * t3 + t30 * t70 + t31 * t69) * t325 + (t4 * t73 + t94 * t47 + t93 * t48) * t326 + t262 * t148 * t132 * t327 + t192 + (t102 * t336 + t337 * t162 + t341 * t260 + t92 * t342 + (-t340 - t347 + t352) * t259) * t166 + (t313 * t342 - t101 * t336 + (t338 * t164 + (-t346 + t347) * qJD(1)) * t164 + t340 * t260 + t339 * t259 + ((t100 * t253 + t252 * t98 + t255 * t88 + t257 * t90 + t337) * t164 + (-t252 * t97 - t253 * t99 - t255 * t87 - t257 * t89 + t338 + t343) * t166 + ((-t298 - t303 - t275 - t295) * t164 + (t299 + t304 + t296 + t297) * t166) * qJD(3) + (t349 + (-t331 - t353) * t166 + t339 + t341) * qJD(1)) * t166) * t164; m(5) * (t164 * t25 + t166 * t24 + (-t164 * t60 + t166 * t59) * qJD(1)) + m(6) * (-qJD(1) * t214 + t14 * t166 + t15 * t164); 0; m(6) * (-qJD(1) * t213 + t164 * t30 + t166 * t31) + m(5) * (t164 * t47 + t166 * t48 + (-t164 * t93 + t166 * t94) * qJD(1)); 0; m(6) * (t214 * t106 + (t14 * t164 - t15 * t166 + (t164 * t49 + t166 * t50) * qJD(1)) * t114) + t184; m(6) * t12; m(6) * (t12 * t7 + t36 * t3 + t213 * t106 + (t164 * t31 - t166 * t30 + (t164 * t70 + t166 * t69) * qJD(1)) * t114) + t192; 0; (t106 * t114 * t262 + t36 * t12) * t325 + t192;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

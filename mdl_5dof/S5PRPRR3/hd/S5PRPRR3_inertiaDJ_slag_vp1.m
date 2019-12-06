% Calculate time derivative of joint inertia matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:39
% EndTime: 2019-12-05 15:46:52
% DurationCPUTime: 7.29s
% Computational Cost: add. (17591->480), mult. (19235->767), div. (0->0), fcn. (18340->10), ass. (0->267)
t207 = qJ(2) + pkin(9);
t202 = sin(t207);
t203 = cos(t207);
t213 = sin(qJ(2));
t215 = cos(qJ(2));
t364 = (-Icges(3,5) * t213 - Icges(4,5) * t202 - Icges(3,6) * t215 - Icges(4,6) * t203) * qJD(2);
t210 = cos(pkin(8));
t347 = t210 ^ 2;
t209 = sin(pkin(8));
t348 = t209 ^ 2;
t363 = t347 + t348;
t353 = qJD(2) * t363;
t362 = t364 * t209;
t361 = t364 * t210;
t214 = cos(qJ(4));
t340 = pkin(4) * t214;
t147 = -pkin(7) * t203 + t202 * t340;
t350 = 2 * m(5);
t349 = 2 * m(6);
t346 = -t203 / 0.2e1;
t345 = t209 / 0.2e1;
t344 = -t210 / 0.2e1;
t343 = t210 / 0.2e1;
t342 = pkin(2) * t213;
t337 = pkin(2) * qJD(2);
t208 = qJ(4) + qJ(5);
t205 = cos(t208);
t310 = t205 * t210;
t204 = sin(t208);
t313 = t204 * t209;
t174 = -t203 * t313 - t310;
t311 = t205 * t209;
t312 = t204 * t210;
t175 = t203 * t311 - t312;
t319 = t202 * t209;
t111 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t319;
t113 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t319;
t115 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t319;
t176 = -t203 * t312 + t311;
t177 = t203 * t310 + t313;
t318 = t202 * t210;
t50 = t111 * t318 + t113 * t176 + t115 * t177;
t336 = t209 * t50;
t304 = t210 * t214;
t212 = sin(qJ(4));
t307 = t209 * t212;
t188 = -t203 * t307 - t304;
t305 = t210 * t212;
t306 = t209 * t214;
t189 = t203 * t306 - t305;
t127 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t319;
t129 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t319;
t131 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t319;
t190 = -t203 * t305 + t306;
t191 = t203 * t304 + t307;
t59 = t127 * t318 + t129 * t190 + t131 * t191;
t335 = t209 * t59;
t112 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t318;
t114 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t318;
t116 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t318;
t49 = t112 * t319 + t114 * t174 + t116 * t175;
t334 = t210 * t49;
t128 = Icges(5,5) * t191 + Icges(5,6) * t190 + Icges(5,3) * t318;
t130 = Icges(5,4) * t191 + Icges(5,2) * t190 + Icges(5,6) * t318;
t132 = Icges(5,1) * t191 + Icges(5,4) * t190 + Icges(5,5) * t318;
t58 = t128 * t319 + t130 * t188 + t132 * t189;
t333 = t210 * t58;
t292 = qJD(4) * t212;
t288 = pkin(4) * t292;
t217 = -qJD(2) * t147 - t203 * t288;
t291 = qJD(4) * t214;
t287 = pkin(4) * t291;
t294 = qJD(2) * t210;
t206 = qJD(4) + qJD(5);
t309 = t206 * t209;
t231 = t202 * t294 - t309;
t308 = t206 * t210;
t284 = t203 * t308;
t139 = t204 * t231 - t205 * t284;
t140 = -t204 * t284 - t205 * t231;
t280 = t203 * t294;
t83 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t280;
t332 = t209 * t287 + t210 * t217 + t83;
t117 = rSges(6,1) * t175 + rSges(6,2) * t174 + rSges(6,3) * t319;
t295 = qJD(2) * t209;
t230 = t202 * t295 + t308;
t285 = t203 * t309;
t137 = t204 * t230 - t205 * t285;
t138 = -t204 * t285 - t205 * t230;
t281 = t203 * t295;
t82 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t281;
t331 = t117 * t280 + t82 * t318;
t326 = Icges(5,4) * t212;
t325 = Icges(5,4) * t214;
t324 = Icges(6,4) * t204;
t323 = Icges(6,4) * t205;
t253 = Icges(6,5) * t205 - Icges(6,6) * t204;
t148 = -Icges(6,3) * t203 + t202 * t253;
t322 = t148 * t203;
t320 = t202 * t206;
t317 = t203 * ((-Icges(6,5) * t204 - Icges(6,6) * t205) * t320 + (Icges(6,3) * t202 + t203 * t253) * qJD(2));
t254 = Icges(5,5) * t214 - Icges(5,6) * t212;
t293 = qJD(4) * t202;
t316 = t203 * ((-Icges(5,5) * t212 - Icges(5,6) * t214) * t293 + (Icges(5,3) * t202 + t203 * t254) * qJD(2));
t156 = -Icges(5,3) * t203 + t202 * t254;
t315 = t203 * t156;
t314 = t203 * t209;
t265 = rSges(6,1) * t205 - rSges(6,2) * t204;
t110 = (-rSges(6,1) * t204 - rSges(6,2) * t205) * t320 + (rSges(6,3) * t202 + t203 * t265) * qJD(2);
t218 = pkin(7) * t202 + t203 * t340;
t135 = qJD(2) * t218 - t202 * t288;
t303 = -t110 - t135;
t124 = -pkin(4) * t305 + t209 * t218;
t302 = t117 + t124;
t118 = rSges(6,1) * t177 + rSges(6,2) * t176 + rSges(6,3) * t318;
t125 = pkin(4) * t307 + t210 * t218;
t301 = t118 + t125;
t151 = -rSges(6,3) * t203 + t202 * t265;
t80 = t203 * t117 + t151 * t319;
t300 = -t147 - t151;
t299 = t363 * pkin(2) * t215;
t298 = t363 * t213 * t337;
t297 = qJD(2) * t202;
t296 = qJD(2) * t203;
t289 = t215 * t337;
t286 = t110 * t319 + t151 * t281 + t203 * t82;
t283 = t212 * t297;
t282 = t214 * t297;
t196 = rSges(4,1) * t202 + rSges(4,2) * t203;
t278 = -t196 - t342;
t197 = t202 * pkin(3) - t203 * pkin(6);
t277 = -t197 - t342;
t276 = t301 * t203;
t275 = t300 * t210;
t269 = pkin(3) * t203 + pkin(6) * t202;
t274 = t269 * t363 + t299;
t273 = -t197 * t353 - t298;
t266 = rSges(5,1) * t214 - rSges(5,2) * t212;
t161 = -rSges(5,3) * t203 + t202 * t266;
t272 = -t161 + t277;
t267 = rSges(4,1) * t203 - rSges(4,2) * t202;
t271 = -t267 * qJD(2) - t289;
t270 = -t269 * qJD(2) - t289;
t199 = rSges(3,1) * t213 + rSges(3,2) * t215;
t252 = -t113 * t204 + t115 * t205;
t55 = -t111 * t203 + t202 * t252;
t251 = -t114 * t204 + t116 * t205;
t56 = -t112 * t203 + t202 * t251;
t264 = t55 * t209 + t56 * t210;
t250 = -t129 * t212 + t131 * t214;
t61 = -t127 * t203 + t202 * t250;
t249 = -t130 * t212 + t132 * t214;
t62 = -t128 * t203 + t202 * t249;
t263 = t61 * t209 + t62 * t210;
t260 = Icges(5,1) * t214 - t326;
t259 = Icges(6,1) * t205 - t324;
t256 = -Icges(5,2) * t212 + t325;
t255 = -Icges(6,2) * t204 + t323;
t133 = rSges(5,1) * t189 + rSges(5,2) * t188 + rSges(5,3) * t319;
t134 = rSges(5,1) * t191 + rSges(5,2) * t190 + rSges(5,3) * t318;
t248 = t133 * t210 - t134 * t209;
t149 = -Icges(6,6) * t203 + t202 * t255;
t150 = -Icges(6,5) * t203 + t202 * t259;
t247 = t149 * t204 - t150 * t205;
t157 = -Icges(5,6) * t203 + t202 * t256;
t158 = -Icges(5,5) * t203 + t202 * t260;
t246 = t157 * t212 - t158 * t214;
t239 = t277 + t300;
t126 = (-rSges(5,1) * t212 - rSges(5,2) * t214) * t293 + (rSges(5,3) * t202 + t203 * t266) * qJD(2);
t238 = -t126 + t270;
t74 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t281;
t237 = t111 * t296 + t202 * t74;
t75 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t280;
t236 = t112 * t296 + t202 * t75;
t143 = -qJD(4) * t189 + t209 * t283;
t144 = qJD(4) * t188 - t209 * t282;
t92 = Icges(5,5) * t144 + Icges(5,6) * t143 + Icges(5,3) * t281;
t235 = t127 * t296 + t202 * t92;
t145 = -qJD(4) * t191 + t210 * t283;
t146 = qJD(4) * t190 - t210 * t282;
t93 = Icges(5,5) * t146 + Icges(5,6) * t145 + Icges(5,3) * t280;
t234 = t128 * t296 + t202 * t93;
t89 = -t196 * t353 - t298;
t232 = t89 * t267;
t107 = (-Icges(6,2) * t205 - t324) * t320 + (Icges(6,6) * t202 + t203 * t255) * qJD(2);
t108 = (-Icges(6,1) * t204 - t323) * t320 + (Icges(6,5) * t202 + t203 * t259) * qJD(2);
t76 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t281;
t78 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t281;
t17 = (qJD(2) * t252 - t74) * t203 + (qJD(2) * t111 + (-t113 * t206 + t78) * t205 + (-t115 * t206 - t76) * t204) * t202;
t77 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t280;
t79 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t280;
t18 = (qJD(2) * t251 - t75) * t203 + (qJD(2) * t112 + (-t114 * t206 + t79) * t205 + (-t116 * t206 - t77) * t204) * t202;
t48 = t111 * t319 + t113 * t174 + t115 * t175;
t51 = t112 * t318 + t114 * t176 + t116 * t177;
t22 = t113 * t137 + t115 * t138 + t174 * t76 + t175 * t78 + t209 * t237;
t23 = t114 * t137 + t116 * t138 + t174 * t77 + t175 * t79 + t209 * t236;
t65 = t148 * t319 + t149 * t174 + t150 * t175;
t6 = -(t107 * t174 + t108 * t175 + t137 * t149 + t138 * t150) * t203 + (t23 * t210 + (t22 - t317) * t209) * t202 + (t65 * t202 + (t334 + (t48 - t322) * t209) * t203) * qJD(2);
t66 = t148 * t318 + t149 * t176 + t150 * t177;
t24 = t113 * t139 + t115 * t140 + t176 * t76 + t177 * t78 + t210 * t237;
t25 = t114 * t139 + t116 * t140 + t176 * t77 + t177 * t79 + t210 * t236;
t7 = -(t107 * t176 + t108 * t177 + t139 * t149 + t140 * t150) * t203 + (t24 * t209 + (t25 - t317) * t210) * t202 + (t66 * t202 + (t336 + (t51 - t322) * t210) * t203) * qJD(2);
t71 = -t202 * t247 - t322;
t229 = -t203 * ((t317 + (t203 * t247 + t264) * qJD(2)) * t203 + (t18 * t210 + t17 * t209 - (qJD(2) * t148 + (-t149 * t206 + t108) * t205 + (-t150 * t206 - t107) * t204) * t203 + t71 * qJD(2)) * t202) + t7 * t318 + t6 * t319 + (-t203 * t65 + (t209 * t48 + t334) * t202) * t281 + (-t203 * t66 + (t210 * t51 + t336) * t202) * t280 + (t202 * t264 - t203 * t71) * t297;
t228 = t270 + t303;
t13 = t209 * t23 - t210 * t22;
t14 = t209 * t25 - t210 * t24;
t227 = t6 * t344 + t7 * t345 + (-t17 * t210 + t18 * t209) * t346 + t13 * t319 / 0.2e1 + t14 * t318 / 0.2e1 + (t209 * t56 - t210 * t55) * t297 / 0.2e1 + (t209 * (t209 * t49 - t210 * t48) + t210 * (t209 * t51 - t210 * t50)) * t296 / 0.2e1;
t153 = t271 * t210;
t152 = t271 * t209;
t136 = t199 * t353;
t123 = (-Icges(5,1) * t212 - t325) * t293 + (Icges(5,5) * t202 + t203 * t260) * qJD(2);
t122 = (-Icges(5,2) * t214 - t326) * t293 + (Icges(5,6) * t202 + t203 * t256) * qJD(2);
t120 = t272 * t210;
t119 = t272 * t209;
t105 = t118 * t297;
t104 = t117 * t318;
t100 = t209 * t217 - t210 * t287;
t99 = rSges(5,1) * t146 + rSges(5,2) * t145 + rSges(5,3) * t280;
t98 = rSges(5,1) * t144 + rSges(5,2) * t143 + rSges(5,3) * t281;
t97 = Icges(5,1) * t146 + Icges(5,4) * t145 + Icges(5,5) * t280;
t96 = Icges(5,1) * t144 + Icges(5,4) * t143 + Icges(5,5) * t281;
t95 = Icges(5,4) * t146 + Icges(5,2) * t145 + Icges(5,6) * t280;
t94 = Icges(5,4) * t144 + Icges(5,2) * t143 + Icges(5,6) * t281;
t91 = t238 * t210;
t90 = t238 * t209;
t88 = -t134 * t203 - t161 * t318;
t87 = t133 * t203 + t161 * t319;
t86 = t239 * t210;
t85 = t239 * t209;
t84 = -t202 * t246 - t315;
t81 = -t118 * t203 - t151 * t318;
t70 = t248 * t202;
t69 = t156 * t318 + t157 * t190 + t158 * t191;
t68 = t156 * t319 + t157 * t188 + t158 * t189;
t67 = -t118 * t319 + t104;
t64 = t228 * t210;
t63 = t228 * t209;
t60 = t128 * t318 + t130 * t190 + t132 * t191;
t57 = t127 * t319 + t129 * t188 + t131 * t189;
t54 = t202 * t275 - t276;
t53 = t124 * t203 + t147 * t319 + t80;
t52 = t133 * t209 + t134 * t210 + t274;
t47 = -t126 * t318 - t203 * t99 + (-t161 * t203 * t210 + t134 * t202) * qJD(2);
t46 = t126 * t319 + t203 * t98 + (-t133 * t202 + t161 * t314) * qJD(2);
t45 = t104 + (t124 * t210 - t209 * t301) * t202;
t44 = t209 * t98 + t210 * t99 + t273;
t43 = -t110 * t318 + t105 + (-t151 * t294 - t83) * t203;
t42 = -t117 * t297 + t286;
t41 = (-t209 * t99 + t210 * t98) * t202 + t248 * t296;
t40 = t209 * t302 + t210 * t301 + t274;
t39 = (-t118 * t296 - t202 * t83) * t209 + t331;
t37 = t332 * t210 + (t100 + t82) * t209 + t273;
t34 = t105 + (qJD(2) * t125 + t210 * t303) * t202 + (qJD(2) * t275 - t332) * t203;
t33 = t135 * t319 + t100 * t203 + (t147 * t314 - t202 * t302) * qJD(2) + t286;
t32 = t130 * t145 + t132 * t146 + t190 * t95 + t191 * t97 + t210 * t234;
t31 = t129 * t145 + t131 * t146 + t190 * t94 + t191 * t96 + t210 * t235;
t30 = t130 * t143 + t132 * t144 + t188 * t95 + t189 * t97 + t209 * t234;
t29 = t129 * t143 + t131 * t144 + t188 * t94 + t189 * t96 + t209 * t235;
t27 = (qJD(2) * t249 - t93) * t203 + (qJD(2) * t128 - t212 * t95 + t214 * t97 + (-t130 * t214 - t132 * t212) * qJD(4)) * t202;
t26 = (qJD(2) * t250 - t92) * t203 + (qJD(2) * t127 - t212 * t94 + t214 * t96 + (-t129 * t214 - t131 * t212) * qJD(4)) * t202;
t19 = (t100 * t202 + t124 * t296) * t210 + (-qJD(2) * t276 - t202 * t332) * t209 + t331;
t16 = t209 * t32 - t210 * t31;
t15 = t209 * t30 - t210 * t29;
t9 = -(t122 * t190 + t123 * t191 + t145 * t157 + t146 * t158) * t203 + (t31 * t209 + (t32 - t316) * t210) * t202 + (t69 * t202 + (t335 + (t60 - t315) * t210) * t203) * qJD(2);
t8 = -(t122 * t188 + t123 * t189 + t143 * t157 + t144 * t158) * t203 + (t30 * t210 + (t29 - t316) * t209) * t202 + (t68 * t202 + (t333 + (t57 - t315) * t209) * t203) * qJD(2);
t1 = [0; -m(3) * t136 + m(4) * t89 + m(5) * t44 + m(6) * t37; (t37 * t40 + t63 * t85 + t64 * t86) * t349 - t210 * t13 - t210 * t15 + (t119 * t90 + t120 * t91 + t44 * t52) * t350 + 0.2e1 * m(4) * (t299 * t89 + (t153 * t278 + t210 * t232) * t210) + 0.2e1 * m(3) * (qJD(2) * t199 - t136) * t363 * (rSges(3,1) * t215 - rSges(3,2) * t213) - t362 * t210 * t347 + (t14 + t16 + t361 * t348 + 0.2e1 * (t152 * t278 + t209 * t232) * m(4) + (-t362 * t209 + t361 * t210) * t210) * t209; 0; m(6) * (t209 * t64 - t210 * t63) + m(5) * (t209 * t91 - t210 * t90) + m(4) * (-t152 * t210 + t153 * t209); 0; m(5) * t41 + m(6) * t19; t9 * t345 + (t27 * t209 - t26 * t210) * t346 + t8 * t344 + (t15 * t345 + t16 * t343) * t202 + m(6) * (t19 * t40 + t33 * t86 + t34 * t85 + t37 * t45 + t53 * t64 + t54 * t63) + m(5) * (t119 * t47 + t120 * t46 + t41 * t52 + t44 * t70 + t87 * t91 + t88 * t90) + (t202 * (t209 * t62 - t210 * t61) / 0.2e1 + ((t60 * t209 - t59 * t210) * t343 + (t58 * t209 - t57 * t210) * t345) * t203) * qJD(2) + t227; m(5) * (t209 * t46 - t210 * t47) + m(6) * (t209 * t33 - t210 * t34); (-t203 * t68 + (t209 * t57 + t333) * t202) * t281 + t8 * t319 + (-t203 * t69 + (t210 * t60 + t335) * t202) * t280 + t9 * t318 + (t19 * t45 + t33 * t53 + t34 * t54) * t349 + (t41 * t70 + t46 * t87 + t47 * t88) * t350 + (t202 * t263 - t84 * t203) * t297 - t203 * ((t316 + (t203 * t246 + t263) * qJD(2)) * t203 + (t27 * t210 + t26 * t209 - (qJD(2) * t156 - t122 * t212 + t123 * t214 - t157 * t291 - t158 * t292) * t203 + t84 * qJD(2)) * t202) + t229; m(6) * t39; m(6) * (t37 * t67 + t39 * t40 + t42 * t86 + t43 * t85 + t63 * t81 + t64 * t80) + t227; m(6) * (t209 * t42 - t210 * t43); m(6) * (t19 * t67 + t33 * t80 + t34 * t81 + t39 * t45 + t42 * t53 + t43 * t54) + t229; (t39 * t67 + t42 * t80 + t43 * t81) * t349 + t229;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

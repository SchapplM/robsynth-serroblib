% Calculate vector of inverse dynamics joint torques for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:46
% EndTime: 2019-12-05 17:28:59
% DurationCPUTime: 9.44s
% Computational Cost: add. (9352->504), mult. (9298->667), div. (0->0), fcn. (8638->10), ass. (0->228)
t181 = pkin(9) + qJ(5);
t176 = sin(t181);
t178 = cos(t181);
t184 = sin(pkin(8));
t186 = cos(pkin(8));
t109 = -rSges(6,3) * t186 + (rSges(6,1) * t178 - rSges(6,2) * t176) * t184;
t168 = -qJD(5) * t186 + qJD(1);
t182 = qJ(1) + pkin(7);
t179 = cos(t182);
t269 = qJD(5) * t184;
t177 = sin(t182);
t289 = t179 * t186;
t112 = t176 * t289 - t177 * t178;
t113 = t176 * t177 + t178 * t289;
t226 = t113 * rSges(6,1) - t112 * rSges(6,2);
t290 = t179 * t184;
t263 = rSges(6,3) * t290;
t67 = -t263 - t226;
t336 = -t109 * t179 * t269 - t168 * t67;
t292 = t177 * t186;
t110 = t176 * t292 + t178 * t179;
t111 = -t179 * t176 + t178 * t292;
t98 = Icges(6,4) * t113;
t62 = Icges(6,2) * t112 - Icges(6,6) * t290 - t98;
t97 = Icges(6,4) * t112;
t64 = Icges(6,1) * t113 + Icges(6,5) * t290 - t97;
t315 = -t110 * t62 - t111 * t64;
t222 = t112 * t62 + t113 * t64;
t58 = Icges(6,5) * t113 - Icges(6,6) * t112 + Icges(6,3) * t290;
t24 = -t186 * t58 + (t176 * t62 + t178 * t64) * t184;
t190 = qJD(1) ^ 2;
t333 = t179 * qJ(3);
t185 = cos(pkin(9));
t286 = t185 * t186;
t183 = sin(pkin(9));
t288 = t183 * t186;
t294 = t177 * t183;
t227 = (t179 * t286 + t294) * rSges(5,1) - (-t177 * t185 + t179 * t288) * rSges(5,2);
t332 = -rSges(5,3) * t290 - t227;
t291 = t179 * t183;
t330 = -rSges(5,2) * (t177 * t288 + t179 * t185) - rSges(5,1) * (-t177 * t286 + t291);
t308 = rSges(3,1) * t179;
t189 = cos(qJ(1));
t320 = pkin(1) * t189;
t234 = -t308 - t320;
t139 = t177 * rSges(3,2) + t234;
t293 = t177 * t184;
t57 = -Icges(6,5) * t111 + Icges(6,6) * t110 - Icges(6,3) * t293;
t299 = Icges(6,4) * t111;
t60 = Icges(6,2) * t110 - Icges(6,6) * t293 - t299;
t96 = Icges(6,4) * t110;
t63 = -Icges(6,1) * t111 - Icges(6,5) * t293 + t96;
t18 = -t112 * t60 + t113 * t63 + t57 * t290;
t197 = t177 * (Icges(6,2) * t111 + t63 + t96) - t179 * (-Icges(6,2) * t113 + t64 - t97);
t198 = t177 * (-Icges(6,1) * t110 - t299 + t60) - t179 * (Icges(6,1) * t112 - t62 + t98);
t323 = -m(5) - m(6);
t267 = qJD(1) * qJD(5);
t126 = (-qJDD(5) * t177 - t179 * t267) * t184;
t322 = t126 / 0.2e1;
t127 = (qJDD(5) * t179 - t177 * t267) * t184;
t321 = t127 / 0.2e1;
t319 = pkin(3) * t186;
t318 = pkin(4) * t183;
t188 = sin(qJ(1));
t317 = t188 * pkin(1);
t167 = -qJDD(5) * t186 + qJDD(1);
t297 = Icges(6,4) * t178;
t107 = -Icges(6,6) * t186 + (-Icges(6,2) * t176 + t297) * t184;
t298 = Icges(6,4) * t176;
t108 = -Icges(6,5) * t186 + (Icges(6,1) * t178 - t298) * t184;
t133 = (-Icges(6,5) * t176 - Icges(6,6) * t178) * t184;
t118 = qJD(5) * t133;
t134 = (-Icges(6,2) * t178 - t298) * t184;
t119 = qJD(5) * t134;
t135 = (-Icges(6,1) * t176 - t297) * t184;
t120 = qJD(5) * t135;
t26 = -t118 * t186 + (-t119 * t176 + t120 * t178 + (-t107 * t178 - t108 * t176) * qJD(5)) * t184;
t106 = -Icges(6,3) * t186 + (Icges(6,5) * t178 - Icges(6,6) * t176) * t184;
t46 = -t106 * t186 + (-t107 * t176 + t108 * t178) * t184;
t316 = t46 * t167 + t26 * t168;
t314 = -t110 * t60 + t111 * t63;
t78 = qJD(1) * t112 + qJD(5) * t111;
t79 = -qJD(1) * t113 + qJD(5) * t110;
t309 = t79 * rSges(6,1) + t78 * rSges(6,2);
t307 = rSges(4,1) * t186;
t306 = rSges(6,3) * t184;
t30 = t106 * t290 - t107 * t112 + t108 * t113;
t305 = t168 * t30;
t23 = -t186 * t57 + (-t176 * t60 + t178 * t63) * t184;
t304 = t23 * t126;
t303 = t24 * t127;
t302 = rSges(4,3) + qJ(3);
t301 = -t111 * rSges(6,1) + t110 * rSges(6,2);
t273 = qJD(1) * t179;
t171 = qJD(3) * t177;
t274 = qJD(1) * t177;
t275 = -pkin(2) * t274 + t171;
t115 = qJ(3) * t273 + t275;
t271 = qJD(4) * t184;
t157 = t179 * t271;
t272 = qJD(1) * t184;
t258 = t177 * t272;
t277 = -qJ(4) * t258 - t274 * t319;
t300 = -t115 - t157 - t277;
t296 = qJ(3) * t177;
t174 = pkin(4) * t185 + pkin(3);
t295 = t174 * t186;
t287 = t184 * qJ(4);
t285 = t189 * t190;
t187 = -pkin(6) - qJ(4);
t284 = qJ(4) + t187;
t282 = t227 * qJD(1);
t281 = t107 - t135;
t280 = t108 + t134;
t145 = pkin(2) * t179 + t296;
t143 = qJD(1) * t145;
t172 = qJD(3) * t179;
t279 = t143 - t172;
t278 = -pkin(4) * t294 - t174 * t289;
t175 = t190 * t317;
t276 = qJDD(3) * t179 + t175;
t136 = (-rSges(6,1) * t176 - rSges(6,2) * t178) * t184;
t121 = qJD(5) * t136;
t270 = qJD(5) * t121;
t268 = -m(4) + t323;
t266 = qJDD(4) * t184;
t265 = pkin(4) * t291;
t262 = qJD(1) * t317;
t261 = rSges(5,3) * t293;
t260 = t109 * t293;
t259 = -t157 - t275;
t257 = t179 * t272;
t256 = t177 * t271;
t254 = -pkin(2) - t307;
t253 = (pkin(3) - t174) * t186;
t251 = -t269 / 0.2e1;
t250 = t269 / 0.2e1;
t249 = -t145 - t320;
t247 = qJ(3) + t318;
t246 = t284 * t184;
t166 = -qJDD(4) * t186 + qJDD(2);
t158 = rSges(4,2) * t290;
t212 = -rSges(4,1) * t289 - rSges(4,3) * t177;
t95 = -t158 - t212;
t245 = t249 - t95;
t244 = t177 * t251;
t243 = t177 * t250;
t242 = t179 * t251;
t241 = t179 * t250;
t223 = t287 + t319;
t137 = t223 * t179;
t240 = -t137 + t249;
t76 = qJD(1) * t110 - qJD(5) * t113;
t77 = -qJD(1) * t111 - qJD(5) * t112;
t239 = rSges(6,1) * t77 + rSges(6,2) * t76;
t34 = Icges(6,5) * t77 + Icges(6,6) * t76 - Icges(6,3) * t258;
t35 = Icges(6,5) * t79 + Icges(6,6) * t78 - Icges(6,3) * t257;
t36 = Icges(6,4) * t77 + Icges(6,2) * t76 - Icges(6,6) * t258;
t37 = Icges(6,4) * t79 + Icges(6,2) * t78 - Icges(6,6) * t257;
t38 = Icges(6,1) * t77 + Icges(6,4) * t76 - Icges(6,5) * t258;
t39 = Icges(6,1) * t79 + Icges(6,4) * t78 - Icges(6,5) * t257;
t238 = -(-t112 * t37 + t113 * t39 + t60 * t76 + t63 * t77 + (t179 * t35 - t274 * t57) * t184) * t177 + t179 * (-t112 * t36 + t113 * t38 - t62 * t76 + t64 * t77 + (t179 * t34 - t274 * t58) * t184);
t237 = -t177 * (t110 * t37 - t111 * t39 + t60 * t78 + t63 * t79 + (-t177 * t35 - t273 * t57) * t184) + t179 * (t110 * t36 - t111 * t38 - t62 * t78 + t64 * t79 + (-t177 * t34 - t273 * t58) * t184);
t236 = -pkin(2) + (-rSges(6,3) + t187) * t184;
t235 = (-rSges(5,3) - qJ(4)) * t184 - pkin(2);
t233 = -t172 + t256;
t10 = -t186 * t34 + (-t176 * t36 + t178 * t38 + (-t176 * t64 + t178 * t62) * qJD(5)) * t184;
t9 = -t186 * t35 + (-t176 * t37 + t178 * t39 + (-t176 * t63 - t178 * t60) * qJD(5)) * t184;
t232 = t10 * t179 - t177 * t9;
t160 = rSges(2,1) * t189 - t188 * rSges(2,2);
t231 = rSges(2,1) * t188 + rSges(2,2) * t189;
t230 = rSges(3,1) * t177 + rSges(3,2) * t179;
t229 = -rSges(4,2) * t184 + t307;
t228 = t330 * qJD(1);
t225 = -t296 - t320;
t224 = -t177 * pkin(2) + t333;
t17 = -t293 * t58 + t315;
t217 = t293 * t57 + t314;
t221 = -t17 * t177 + t179 * t217;
t40 = -rSges(6,3) * t258 + t239;
t41 = -rSges(6,3) * t257 + t309;
t220 = -t177 * t40 - t179 * t41;
t219 = t177 * (Icges(6,5) * t110 + Icges(6,6) * t111) - t179 * (-Icges(6,5) * t112 - Icges(6,6) * t113);
t218 = -t184 * t187 + t295;
t74 = (t246 + t319) * t179 + t278;
t216 = t240 + t74;
t215 = t240 + t332;
t214 = rSges(6,3) * t293 - t301;
t213 = -qJD(1) * t224 - t171 + t262;
t209 = -qJD(1) * t318 - t271;
t205 = (t17 * t179 + t177 * t217) * t184;
t19 = t58 * t290 + t222;
t204 = (-t177 * t18 + t179 * t19) * t184;
t203 = qJD(1) * t223;
t202 = qJD(1) * t137 + t143 + t233;
t200 = t230 + t317;
t199 = t177 * t203 - t157 + t213;
t195 = rSges(4,3) * t179 - t229 * t177;
t194 = qJDD(1) * t224 + qJDD(3) * t177 + (-qJDD(1) * t188 - t285) * pkin(1) + (-t279 + t172) * qJD(1);
t193 = t184 * (-t179 * t301 + (t67 + t263) * t177);
t192 = -qJDD(1) * t223 * t177 + t179 * t266 + qJD(1) * (-t179 * t203 - t256) + t194;
t191 = -t261 - t330;
t164 = rSges(3,2) * t274;
t152 = rSges(4,2) * t257;
t148 = t187 * t257;
t124 = t230 * qJD(1) + t262;
t87 = -rSges(6,1) * t112 - rSges(6,2) * t113;
t86 = rSges(6,1) * t110 + rSges(6,2) * t111;
t69 = qJD(1) * t245 + t172;
t68 = -qJD(1) * t195 + t213;
t43 = qJD(1) * t215 - t233;
t32 = t245 * qJDD(1) + (-rSges(4,3) * t273 - t115 + (qJD(1) * t229 - qJD(3)) * t177) * qJD(1) + t276;
t31 = qJDD(1) * t195 + (qJD(1) * t212 + t152) * qJD(1) + t194;
t29 = -t106 * t293 + t107 * t110 - t108 * t111;
t28 = t29 * t168;
t27 = -qJD(4) * t186 + qJD(5) * t193 + qJD(2);
t22 = qJD(1) * t216 - t233 - t336;
t21 = t168 * t214 - qJD(1) * (t265 + (t253 + t246) * t177) - qJD(5) * t260 + t199;
t15 = -t177 * t266 + t215 * qJDD(1) + (-t157 + (rSges(5,3) * t272 - qJD(3)) * t177 + t228 + t300) * qJD(1) + t276;
t14 = qJDD(1) * t191 + ((-rSges(5,3) * t273 - qJD(4) * t177) * t184 - t282) * qJD(1) + t192;
t13 = t107 * t78 + t108 * t79 + t110 * t119 - t111 * t120 + (-t106 * t273 - t118 * t177) * t184;
t12 = t107 * t76 + t108 * t77 - t112 * t119 + t113 * t120 + (-t106 * t274 + t118 * t179) * t184;
t11 = -t126 * t67 + t127 * t214 + t220 * t269 + t166;
t8 = t127 * t109 + t167 * t67 - t168 * t40 + (-qJDD(4) * t177 + t179 * t270) * t184 + t216 * qJDD(1) + (t209 * t179 + (qJD(1) * t218 - qJD(3)) * t177 + t277 + t300) * qJD(1) + t276;
t7 = qJDD(1) * t265 + t167 * t301 + t168 * t41 - t126 * t109 + (qJDD(1) * t253 + (-t167 * rSges(6,3) + qJDD(1) * t284 + t270) * t184) * t177 + (t148 + t209 * t177 + (t253 + t287) * t273) * qJD(1) + t192;
t6 = qJD(5) * t204 + t305;
t5 = qJD(5) * t205 + t28;
t1 = [t304 / 0.2e1 + t303 / 0.2e1 + t316 + (t28 + (t315 * t179 + (-t19 + t217 + t222) * t177) * t269) * t242 + t30 * t321 + t29 * t322 + (-t305 + (-(-t18 - t315) * t177 + (-t222 - t314) * t179 + (-t58 * t177 ^ 2 + (-t177 * t57 - t179 * t58) * t179) * t184 + t221) * t269 + t6) * t243 + (-t124 * t164 + (qJD(1) * t139 * t200 - t124 * t234) * qJD(1) + (t190 * t230 - g(2) + t175) * t139 + (pkin(1) * t285 + (qJD(1) * t308 - t164) * qJD(1) + g(3)) * t200) * m(3) + (g(2) * t160 + g(3) * t231) * m(2) + (t13 + t9) * t244 + (t10 + t12 + t5) * t241 + (-(t22 + (-t74 + t320) * qJD(1) + t202 + t336) * t21 + t22 * (-t239 + t259) - t21 * (t148 - t233 + t309) + ((t188 * t22 + t189 * t21) * pkin(1) + (-t22 * t247 - t21 * (-pkin(2) - t295 - t306)) * t179 + (t22 * (t218 + t306) + t21 * t247) * t177) * qJD(1) + (t7 - g(3)) * (-t317 + t247 * t179 + (t236 - t295) * t177 + t301) + (t8 - g(2)) * (t179 * t236 + t225 - t226 + t278)) * m(6) + ((-g(3) + t14) * (-t317 + (rSges(5,1) * t183 + rSges(5,2) * t185 + qJ(3)) * t179 + ((-rSges(5,1) * t185 + rSges(5,2) * t183 - pkin(3)) * t186 + t235) * t177) + (t228 + t259 - t277 + (t261 + t317 - t333) * qJD(1)) * t43 + (-t43 - t202 + t233 + t282 + (-t320 - t225 + (rSges(5,3) * t184 + pkin(2) + t223) * t179 + t332) * qJD(1)) * (-qJD(1) * t191 + t199) + (-g(2) + t15) * (-t227 + t179 * (t235 - t319) + t225)) * m(5) + (-(t69 + (t95 + t320) * qJD(1) + t279) * t68 - t69 * t275 - t68 * (t152 + t172) + ((t188 * t69 + t189 * t68) * pkin(1) + (-t254 * t68 - t302 * t69) * t179 + (t229 * t69 + t302 * t68) * t177) * qJD(1) + (t32 - g(2)) * (-t302 * t177 + t254 * t179 + t158 - t320) + (t31 - g(3)) * (-t317 + t302 * t179 + (-pkin(2) - t229) * t177)) * m(4) + (m(3) * (t139 ^ 2 + t200 ^ 2) + m(2) * (t160 ^ 2 + t231 ^ 2) + Icges(2,3) + Icges(3,3) + (Icges(5,3) + Icges(4,2)) * t186 ^ 2 + ((Icges(5,1) * t185 ^ 2 + (-0.2e1 * Icges(5,4) * t185 + Icges(5,2) * t183) * t183 + Icges(4,1)) * t184 + 0.2e1 * (-Icges(5,5) * t185 + Icges(5,6) * t183 + Icges(4,4)) * t186) * t184) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t166 + m(6) * t11 + (-m(3) + t268) * g(1); t268 * (g(2) * t179 + g(3) * t177) + m(4) * (t177 * t31 + t179 * t32) + m(5) * (t14 * t177 + t15 * t179) + m(6) * (t177 * t7 + t8 * t179); t323 * (-g(1) * t186 + (-g(2) * t177 + g(3) * t179) * t184) + m(5) * (t14 * t290 - t15 * t293 - t166 * t186) + m(6) * (-t11 * t186 + t290 * t7 - t293 * t8); -t186 * (t232 * t269 + t303 + t304 + t316) / 0.2e1 + t167 * (-t186 * t46 + (-t177 * t23 + t179 * t24) * t184) / 0.2e1 + t168 * (-t186 * t26 + ((-t24 * t177 - t179 * t23) * qJD(1) + t232) * t184) / 0.2e1 - (-t126 * t217 + t127 * t17 + t13 * t168 + t167 * t29 + t237 * t269) * t293 / 0.2e1 + (-t186 * t29 + t205) * t322 + (-t13 * t186 + (qJD(1) * t221 + t237) * t184) * t244 + (t12 * t168 + t126 * t18 + t127 * t19 + t167 * t30 + t238 * t269) * t290 / 0.2e1 + (-t186 * t30 + t204) * t321 + (-t12 * t186 + ((-t19 * t177 - t18 * t179) * qJD(1) + t238) * t184) * t241 - t168 * (-t186 * t133 * t168 + ((-t176 * t280 - t178 * t281) * t168 + ((t176 * t197 + t178 * t198) * t184 + t219 * t186) * qJD(5)) * t184) / 0.2e1 + ((t110 * t280 + t111 * t281 - t133 * t293) * t168 + (-t110 * t197 - t111 * t198 + t219 * t293) * t269) * t243 + ((-t112 * t280 - t113 * t281 + t133 * t290) * t168 + (t112 * t197 + t113 * t198 - t219 * t290) * t269) * t242 - (t177 * t6 + t179 * t5) * t272 / 0.2e1 + (t11 * t193 + t27 * ((-t177 * t214 + t179 * t67) * qJD(1) + t220) * t184 + t8 * (t109 * t290 - t186 * t67) + t22 * (t186 * t40 + (-t109 * t274 + t121 * t179) * t184) + t7 * (t186 * t214 + t260) - t21 * (-t186 * t41 + (t109 * t273 + t121 * t177) * t184) - (-t21 * t86 - t22 * t87) * t168 - (t27 * (-t177 * t87 - t179 * t86) + (-t177 * t21 + t179 * t22) * t136) * t269 - g(1) * t136 - g(2) * t86 - g(3) * t87) * m(6);];
tau = t1;

% Calculate vector of inverse dynamics joint torques for
% S5RPPPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:48
% DurationCPUTime: 7.59s
% Computational Cost: add. (6235->509), mult. (16358->640), div. (0->0), fcn. (17086->8), ass. (0->225)
t233 = cos(pkin(8));
t234 = cos(pkin(7));
t326 = t233 * t234;
t186 = qJD(5) * t326 + qJD(1);
t231 = sin(pkin(8));
t236 = sin(qJ(1));
t325 = t233 * t236;
t232 = sin(pkin(7));
t238 = cos(qJ(1));
t327 = t232 * t238;
t170 = t231 * t327 + t325;
t235 = sin(qJ(5));
t237 = cos(qJ(5));
t321 = t234 * t238;
t126 = -t170 * t235 + t237 * t321;
t127 = t170 * t237 + t235 * t321;
t324 = t233 * t238;
t293 = t232 * t324;
t169 = t231 * t236 - t293;
t64 = Icges(6,5) * t127 + Icges(6,6) * t126 + Icges(6,3) * t169;
t335 = Icges(6,4) * t127;
t67 = Icges(6,2) * t126 + Icges(6,6) * t169 + t335;
t120 = Icges(6,4) * t126;
t70 = Icges(6,1) * t127 + Icges(6,5) * t169 + t120;
t14 = t126 * t67 + t127 * t70 + t169 * t64;
t328 = t232 * t236;
t172 = -t231 * t328 + t324;
t322 = t234 * t237;
t130 = t172 * t235 + t236 * t322;
t323 = t234 * t236;
t131 = -t172 * t237 + t235 * t323;
t171 = t231 * t238 + t232 * t325;
t66 = Icges(6,5) * t131 + Icges(6,6) * t130 - Icges(6,3) * t171;
t334 = Icges(6,4) * t131;
t68 = -Icges(6,2) * t130 + Icges(6,6) * t171 - t334;
t121 = Icges(6,4) * t130;
t71 = -Icges(6,1) * t131 + Icges(6,5) * t171 - t121;
t15 = -t126 * t68 - t127 * t71 + t169 * t66;
t257 = t14 * t169 - t15 * t171;
t167 = t231 * t234 * t235 + t232 * t237;
t168 = t231 * t322 - t232 * t235;
t96 = -Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t326;
t333 = Icges(6,4) * t168;
t97 = Icges(6,2) * t167 + Icges(6,6) * t326 - t333;
t155 = Icges(6,4) * t167;
t98 = -Icges(6,1) * t168 + Icges(6,5) * t326 + t155;
t31 = t126 * t97 + t127 * t98 + t169 * t96;
t5 = qJD(5) * t257 + t31 * t186;
t307 = qJD(1) * t236;
t150 = -qJD(1) * t293 + t231 * t307;
t151 = t170 * qJD(1);
t304 = qJD(4) * t234;
t305 = qJD(3) * t232;
t363 = (t304 + t305) * t236;
t384 = -pkin(4) * t151 - pkin(6) * t150 - t363;
t382 = (-pkin(3) - qJ(2)) * t236;
t178 = pkin(3) * t238 - qJ(4) * t323;
t202 = t238 * t304;
t329 = qJ(3) * t232;
t258 = pkin(2) * t234 + t329;
t173 = t258 * t236;
t203 = t238 * t305;
t225 = t238 * qJ(2);
t187 = pkin(1) * t236 - t225;
t222 = qJD(2) * t236;
t308 = -qJD(1) * t187 + t222;
t271 = -qJD(1) * t173 + t203 + t308;
t306 = qJD(1) * t238;
t309 = qJ(2) * t306 + t222;
t287 = t203 + t309;
t313 = pkin(3) * t306 + t202;
t381 = -qJD(1) * t178 - t202 - t271 + t287 + t313;
t301 = qJD(5) * t171;
t75 = rSges(6,1) * t131 + rSges(6,2) * t130 - t171 * rSges(6,3);
t99 = -rSges(6,1) * t168 + rSges(6,2) * t167 + rSges(6,3) * t326;
t380 = t186 * t75 + t99 * t301;
t16 = t130 * t67 + t131 * t70 - t171 * t64;
t17 = -t130 * t68 - t131 * t71 - t171 * t66;
t256 = t16 * t169 - t17 * t171;
t21 = -t167 * t68 + t168 * t71 + t326 * t66;
t268 = t178 + t225;
t318 = t172 * rSges(5,1) - t171 * rSges(5,2);
t371 = t268 + t318;
t370 = -rSges(5,3) * t323 + t318;
t118 = t170 * pkin(4) + pkin(6) * t169;
t177 = t236 * pkin(3) + qJ(4) * t321;
t174 = pkin(2) * t321 + qJ(3) * t327;
t189 = t238 * pkin(1) + t236 * qJ(2);
t270 = t174 + t189;
t358 = t270 + t177;
t368 = t118 + t358;
t360 = t172 * pkin(4) + pkin(6) * t171;
t100 = t170 * rSges(5,1) - t169 * rSges(5,2) + rSges(5,3) * t321;
t367 = t358 + t100;
t32 = t130 * t97 + t131 * t98 - t171 * t96;
t356 = t234 ^ 2;
t144 = t236 * rSges(4,1) - rSges(4,2) * t321 + rSges(4,3) * t327;
t362 = t144 + t270;
t277 = -pkin(1) - t329;
t347 = -pkin(2) - qJ(4);
t359 = (-rSges(5,3) + t347) * t234 + t277;
t143 = rSges(3,1) * t321 - rSges(3,2) * t327 + t236 * rSges(3,3);
t223 = qJD(2) * t238;
t242 = -t223 + t363;
t357 = t169 * (-Icges(6,1) * t126 + t335 + t67) - t171 * (-Icges(6,1) * t130 + t334 - t68);
t354 = -m(5) - m(6);
t152 = t171 * qJD(1);
t110 = qJD(5) * t152 + qJDD(5) * t169;
t353 = t110 / 0.2e1;
t111 = qJD(5) * t150 - qJDD(5) * t171;
t352 = t111 / 0.2e1;
t351 = t236 / 0.2e1;
t350 = -t238 / 0.2e1;
t348 = g(1) * t236;
t185 = qJDD(5) * t326 + qJDD(1);
t148 = t167 * qJD(5);
t149 = t168 * qJD(5);
t106 = Icges(6,5) * t148 + Icges(6,6) * t149;
t107 = Icges(6,4) * t148 + Icges(6,2) * t149;
t108 = Icges(6,1) * t148 + Icges(6,4) * t149;
t19 = t106 * t326 + t107 * t167 - t108 * t168 + t148 * t98 + t149 * t97;
t35 = t167 * t97 - t168 * t98 + t326 * t96;
t345 = t35 * t185 + t19 * t186;
t20 = t167 * t67 - t168 * t70 + t326 * t64;
t339 = t20 * t110;
t338 = t21 * t111;
t337 = Icges(6,2) * t168 + t155 + t98;
t336 = Icges(6,1) * t167 + t333 - t97;
t166 = qJD(1) * t189 - t223;
t283 = t236 * t305;
t320 = -t258 * t306 - t166 - t283;
t153 = t172 * qJD(1);
t319 = t153 * rSges(5,1) - t152 * rSges(5,2);
t133 = t143 + t189;
t317 = -t173 - t187;
t311 = rSges(3,2) * t328 + t238 * rSges(3,3);
t142 = rSges(3,1) * t323 - t311;
t316 = -t187 - t142;
t286 = t232 * t307;
t315 = rSges(3,2) * t286 + rSges(3,3) * t306;
t285 = t234 * t307;
t314 = rSges(4,1) * t306 + rSges(4,2) * t285;
t312 = t203 + t222;
t299 = qJD(1) * qJD(2);
t310 = qJDD(2) * t236 + t238 * t299;
t303 = qJD(4) * t236;
t302 = qJD(5) * t169;
t300 = -m(4) + t354;
t298 = qJDD(3) * t232;
t297 = qJDD(4) * t234;
t84 = -qJD(5) * t127 - t153 * t235 - t237 * t285;
t85 = qJD(5) * t126 + t153 * t237 - t235 * t285;
t43 = t85 * rSges(6,1) + t84 * rSges(6,2) + t152 * rSges(6,3);
t73 = t127 * rSges(6,1) + t126 * rSges(6,2) + t169 * rSges(6,3);
t291 = t178 + t317;
t275 = rSges(4,1) * t238 - rSges(4,3) * t328;
t145 = rSges(4,2) * t323 + t275;
t290 = t145 + t317;
t289 = t238 * t298 + t310;
t288 = t202 + t312;
t284 = t234 * t306;
t282 = -rSges(3,1) * t234 - pkin(1);
t280 = t302 / 0.2e1;
t279 = -t301 / 0.2e1;
t278 = t301 / 0.2e1;
t276 = t153 * pkin(4) + pkin(6) * t152;
t183 = -qJDD(3) * t234 + qJDD(4) * t232;
t274 = t291 + t370;
t273 = t291 + t360;
t272 = t238 * t297 + t289;
t267 = (rSges(4,2) - pkin(2)) * t234 - pkin(1);
t265 = g(1) * t238 + g(2) * t236;
t82 = -qJD(5) * t131 - t151 * t235 + t237 * t284;
t83 = qJD(5) * t130 + t151 * t237 + t235 * t284;
t36 = Icges(6,5) * t83 + Icges(6,6) * t82 + Icges(6,3) * t150;
t37 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t152;
t38 = Icges(6,4) * t83 + Icges(6,2) * t82 + Icges(6,6) * t150;
t39 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t152;
t40 = Icges(6,1) * t83 + Icges(6,4) * t82 + Icges(6,5) * t150;
t41 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t152;
t264 = (t130 * t39 + t131 * t41 + t150 * t64 - t171 * t37 + t67 * t82 + t70 * t83) * t169 - t171 * (t130 * t38 + t131 * t40 + t150 * t66 - t171 * t36 - t68 * t82 - t71 * t83);
t263 = t169 * (t126 * t39 + t127 * t41 + t152 * t64 + t169 * t37 + t67 * t84 + t70 * t85) - t171 * (t126 * t38 + t127 * t40 + t152 * t66 + t169 * t36 - t68 * t84 - t71 * t85);
t7 = t148 * t70 + t149 * t67 + t167 * t39 - t168 * t41 + t326 * t37;
t8 = -t148 * t71 - t149 * t68 + t167 * t38 - t168 * t40 + t326 * t36;
t262 = t7 * t169 - t8 * t171;
t261 = t223 - t283;
t190 = rSges(2,1) * t238 - rSges(2,2) * t236;
t188 = rSges(2,1) * t236 + rSges(2,2) * t238;
t260 = -t151 * rSges(5,1) + t150 * rSges(5,2);
t42 = rSges(6,1) * t83 + rSges(6,2) * t82 + rSges(6,3) * t150;
t255 = t169 * t42 + t171 * t43;
t254 = t169 * t75 + t171 * t73;
t253 = t169 * (Icges(6,5) * t126 - Icges(6,6) * t127) - t171 * (Icges(6,5) * t130 - Icges(6,6) * t131);
t250 = -pkin(1) - t258;
t249 = -t283 + t320;
t248 = -qJDD(2) * t238 + qJD(1) * (-pkin(1) * t307 + t309) + qJDD(1) * t189 + t236 * t299;
t243 = (-Icges(6,2) * t127 + t120 + t70) * t169 - (-Icges(6,2) * t131 + t121 - t71) * t171;
t241 = qJDD(1) * t174 + t236 * t298 + t248 + (-t258 * t307 + 0.2e1 * t203) * qJD(1);
t240 = qJDD(1) * t177 + t236 * t297 + t241 + (-qJ(4) * t285 + t202 + t313) * qJD(1);
t137 = qJD(1) * t177 + t234 * t303;
t117 = rSges(6,1) * t167 + rSges(6,2) * t168;
t114 = Icges(6,5) * t167 + Icges(6,6) * t168;
t113 = qJD(1) * t133 - t223;
t112 = qJD(1) * t316 + t222;
t109 = rSges(6,1) * t148 + rSges(6,2) * t149;
t95 = rSges(6,1) * t130 - rSges(6,2) * t131;
t94 = rSges(6,1) * t126 - rSges(6,2) * t127;
t87 = qJD(1) * t362 - t261;
t86 = qJD(1) * t290 + t312;
t63 = qJDD(1) * t143 + qJD(1) * (-rSges(3,1) * t285 + t315) + t248;
t62 = t316 * qJDD(1) + (-qJD(1) * t143 - t166) * qJD(1) + t310;
t57 = qJD(1) * t274 + t288;
t45 = qJDD(1) * t144 + qJD(1) * (-rSges(4,3) * t286 + t314) + t241;
t44 = t290 * qJDD(1) + (-qJD(1) * t144 + t249) * qJD(1) + t289;
t33 = -qJD(3) * t234 + qJD(4) * t232 + qJD(5) * t254;
t29 = qJD(1) * t368 + t186 * t73 - t99 * t302 + t242;
t28 = qJD(1) * t273 + t288 - t380;
t27 = qJD(1) * (-rSges(5,3) * t285 + t319) + t240 + qJDD(1) * t100;
t26 = t274 * qJDD(1) + (-t137 + (-rSges(5,3) * t306 - t303) * t234 + t249 + t260) * qJD(1) + t272;
t13 = t106 * t169 + t107 * t126 + t108 * t127 + t152 * t96 + t84 * t97 + t85 * t98;
t12 = -t106 * t171 + t107 * t130 + t108 * t131 + t150 * t96 + t82 * t97 + t83 * t98;
t11 = qJD(5) * t255 + t110 * t75 - t111 * t73 + t183;
t10 = qJD(1) * t276 + qJDD(1) * t118 - t109 * t302 - t110 * t99 + t185 * t73 + t186 * t43 + t240;
t9 = -t109 * t301 + t111 * t99 - t185 * t75 - t186 * t42 + t273 * qJDD(1) + (-t137 + t320 + t384) * qJD(1) + t272;
t1 = [t345 - m(2) * (-g(1) * t188 + g(2) * t190) + t32 * t352 + t31 * t353 + t338 / 0.2e1 + t339 / 0.2e1 + t5 * t278 + (t13 + t7) * t280 + (t12 + t8 + t5) * t279 + ((-g(2) + t10) * (t73 + t368) + (-g(1) + t9) * (-t75 + t268 + t360) + (t9 * t236 - t348) * t250 + (t223 - t42 + ((t234 * t347 + t277) * t238 + t382) * qJD(1) + t384) * t28 + (t28 + t276 + t43 + (-t360 + (-qJ(4) * t234 + t250) * t236) * qJD(1) + t380 + t381) * t29) * m(6) + (-g(1) * t371 - ((-rSges(5,3) - pkin(2)) * t234 + t277) * t348 + (t27 - g(2)) * t367 + ((-rSges(5,3) * t234 + t250) * t236 + t371) * t26 + (t260 + (t238 * t359 + t382) * qJD(1) - t242) * t57 + (t319 + t57 + (t236 * t359 - t370) * qJD(1) + t381) * (qJD(1) * t367 + t242)) * m(5) + (-(qJD(1) * t145 + t271 - t86) * t87 + t86 * t261 + t87 * (t287 + t314) + (t86 * ((-rSges(4,3) - qJ(3)) * t232 + t267) * t238 + (t86 * (-rSges(4,1) - qJ(2)) + t87 * (-rSges(4,3) * t232 + t250)) * t236) * qJD(1) + (-g(2) + t45) * t362 + (-g(1) + t44) * (t225 + (t267 - t329) * t236 + t275)) * m(4) + (-(-qJD(1) * t142 - t112 + t308) * t113 + t112 * t223 + t113 * (t309 + t315) + (t112 * (rSges(3,2) * t232 + t282) * t238 + (t112 * (-rSges(3,3) - qJ(2)) + t113 * t282) * t236) * qJD(1) + (-g(2) + t63) * t133 + (-g(1) + t62) * (t236 * t282 + t225 + t311)) * m(3) + (m(2) * (t188 ^ 2 + t190 ^ 2) + Icges(2,3) + ((Icges(5,3) + Icges(3,1) + Icges(4,2)) * t232 + 0.2e1 * (-Icges(5,5) * t231 - Icges(5,6) * t233 + Icges(3,4) + Icges(4,6)) * t234) * t232 + (t233 ^ 2 * Icges(5,2) + (Icges(5,1) * t231 + 0.2e1 * Icges(5,4) * t233) * t231 + Icges(3,2) + Icges(4,3)) * t356) * qJDD(1); (-m(3) + t300) * (-g(2) * t238 + t348) + 0.2e1 * (t10 * t350 + t351 * t9) * m(6) + 0.2e1 * (t26 * t351 + t27 * t350) * m(5) + 0.2e1 * (t350 * t45 + t351 * t44) * m(4) + 0.2e1 * (t350 * t63 + t351 * t62) * m(3); t300 * (-g(3) * t234 + t232 * t265) + m(4) * (qJDD(3) * t356 + t327 * t44 + t328 * t45) + m(5) * (-t183 * t234 + t26 * t327 + t27 * t328) + m(6) * (t10 * t328 - t11 * t234 + t327 * t9); t354 * (g(3) * t232 + t234 * t265) + m(5) * (t183 * t232 + t26 * t321 + t27 * t323) + m(6) * (t10 * t323 + t11 * t232 + t321 * t9); t152 * t5 / 0.2e1 + t169 * (t263 * qJD(5) + t110 * t14 + t111 * t15 + t13 * t186 + t185 * t31) / 0.2e1 + (t31 * t326 + t257) * t353 + (t13 * t326 + t14 * t152 + t15 * t150 + t263) * t280 + t150 * (qJD(5) * t256 + t186 * t32) / 0.2e1 - t171 * (t264 * qJD(5) + t110 * t16 + t111 * t17 + t12 * t186 + t185 * t32) / 0.2e1 + (t32 * t326 + t256) * t352 + (t12 * t326 + t150 * t17 + t152 * t16 + t264) * t279 + (t262 * qJD(5) + t338 + t339 + t345) * t326 / 0.2e1 + t185 * (t169 * t20 - t171 * t21 + t35 * t326) / 0.2e1 + t186 * (t21 * t150 + t20 * t152 + t19 * t326 + t262) / 0.2e1 - ((t114 * t169 + t126 * t337 + t127 * t336) * t186 + (t126 * t243 - t127 * t357 + t169 * t253) * qJD(5)) * t302 / 0.2e1 + ((-t114 * t171 + t130 * t337 + t131 * t336) * t186 + (t130 * t243 - t131 * t357 - t171 * t253) * qJD(5)) * t278 - t186 * ((t114 * t326 + t337 * t167 - t336 * t168) * t186 + (t243 * t167 + t168 * t357 + t253 * t326) * qJD(5)) / 0.2e1 + (t9 * (-t171 * t99 - t326 * t75) + t28 * (-t109 * t171 + t150 * t99 - t326 * t42) + t10 * (-t169 * t99 + t326 * t73) + t29 * (-t109 * t169 - t152 * t99 + t326 * t43) + t11 * t254 + t33 * (-t150 * t73 + t152 * t75 + t255) - (-t28 * t95 + t29 * t94) * t186 - (t33 * (t169 * t95 + t171 * t94) + (-t169 * t29 - t171 * t28) * t117) * qJD(5) - g(1) * t94 - g(2) * t95 - g(3) * t117) * m(6);];
tau = t1;

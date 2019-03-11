% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:10
% EndTime: 2019-03-08 19:39:37
% DurationCPUTime: 18.02s
% Computational Cost: add. (8280->627), mult. (20015->861), div. (0->0), fcn. (16149->18), ass. (0->287)
t240 = cos(pkin(12));
t226 = pkin(5) * t240 + pkin(4);
t235 = pkin(12) + qJ(6);
t229 = sin(t235);
t231 = cos(t235);
t237 = sin(pkin(12));
t273 = -t240 * mrSges(6,1) + t237 * mrSges(6,2);
t375 = m(6) * pkin(4) + m(7) * t226 + t231 * mrSges(7,1) - t229 * mrSges(7,2) + mrSges(5,1) - t273;
t339 = pkin(9) + qJ(5);
t374 = -m(6) * qJ(5) - m(7) * t339 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t238 = sin(pkin(11));
t246 = sin(qJ(4));
t241 = cos(pkin(11));
t344 = cos(qJ(4));
t292 = t344 * t241;
t260 = -t246 * t238 + t292;
t193 = t260 * qJD(4);
t201 = t238 * t344 + t246 * t241;
t194 = t201 * qJD(4);
t247 = sin(qJ(2));
t239 = sin(pkin(6));
t309 = qJD(1) * t239;
t289 = t247 * t309;
t423 = pkin(4) * t194 - qJ(5) * t193 - qJD(5) * t201 - t289;
t249 = cos(qJ(2));
t315 = t239 * t249;
t253 = t260 * t315;
t340 = pkin(8) + qJ(3);
t209 = t340 * t238;
t211 = t340 * t241;
t383 = -t344 * t209 - t246 * t211;
t386 = -qJD(1) * t253 + t260 * qJD(3) + qJD(4) * t383;
t189 = t260 * qJD(2);
t245 = sin(qJ(6));
t248 = cos(qJ(6));
t200 = t237 * t248 + t240 * t245;
t112 = t200 * t189;
t192 = t200 * qJD(6);
t422 = t112 - t192;
t263 = t237 * t245 - t240 * t248;
t113 = t263 * t189;
t191 = t263 * qJD(6);
t421 = t113 - t191;
t236 = pkin(11) + qJ(4);
t230 = sin(t236);
t232 = cos(t236);
t420 = t230 * t374 - t232 * t375;
t182 = qJD(6) - t189;
t190 = t201 * qJD(2);
t163 = qJD(4) * t237 + t190 * t240;
t276 = t240 * qJD(4) - t190 * t237;
t89 = t163 * t248 + t245 * t276;
t343 = Ifges(7,4) * t89;
t411 = -t163 * t245 + t248 * t276;
t38 = Ifges(7,2) * t411 + Ifges(7,6) * t182 + t343;
t364 = t38 / 0.2e1;
t86 = Ifges(7,4) * t411;
t39 = Ifges(7,1) * t89 + Ifges(7,5) * t182 + t86;
t363 = t39 / 0.2e1;
t395 = -t237 * t386 + t240 * t423;
t394 = t237 * t423 + t240 * t386;
t319 = t193 * t240;
t418 = pkin(5) * t194 - pkin(9) * t319 + t395;
t320 = t193 * t237;
t417 = pkin(9) * t320 - t394;
t414 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t413 = t238 ^ 2 + t241 ^ 2;
t272 = mrSges(6,1) * t237 + mrSges(6,2) * t240;
t412 = -m(7) * pkin(5) * t237 - t229 * mrSges(7,1) - t231 * mrSges(7,2) - mrSges(5,3) - t272;
t143 = qJD(2) * t193 + qJDD(2) * t201;
t114 = qJDD(4) * t240 - t143 * t237;
t115 = qJDD(4) * t237 + t143 * t240;
t33 = qJD(6) * t411 + t114 * t245 + t115 * t248;
t366 = t33 / 0.2e1;
t34 = -qJD(6) * t89 + t114 * t248 - t115 * t245;
t365 = t34 / 0.2e1;
t205 = qJD(2) * qJ(3) + t289;
t242 = cos(pkin(6));
t308 = qJD(1) * t242;
t220 = t241 * t308;
t333 = pkin(8) * qJD(2);
t151 = t220 + (-t205 - t333) * t238;
t167 = t241 * t205 + t238 * t308;
t152 = t241 * t333 + t167;
t84 = t246 * t151 + t152 * t344;
t80 = qJD(4) * qJ(5) + t84;
t227 = pkin(3) * t241 + pkin(2);
t267 = -t249 * t309 + qJD(3);
t178 = -qJD(2) * t227 + t267;
t96 = -pkin(4) * t189 - qJ(5) * t190 + t178;
t40 = -t237 * t80 + t240 * t96;
t23 = -pkin(5) * t189 - pkin(9) * t163 + t40;
t41 = t237 * t96 + t240 * t80;
t29 = pkin(9) * t276 + t41;
t8 = t23 * t248 - t245 * t29;
t410 = t8 * mrSges(7,1);
t9 = t23 * t245 + t248 * t29;
t409 = t9 * mrSges(7,2);
t357 = t114 / 0.2e1;
t356 = t115 / 0.2e1;
t300 = qJDD(2) * t238;
t144 = qJD(2) * t194 - qJDD(2) * t292 + t246 * t300;
t138 = qJDD(6) + t144;
t355 = t138 / 0.2e1;
t354 = t144 / 0.2e1;
t317 = t201 * t240;
t139 = -pkin(4) * t260 - qJ(5) * t201 - t227;
t158 = -t246 * t209 + t211 * t344;
t77 = t240 * t139 - t158 * t237;
t56 = -pkin(5) * t260 - pkin(9) * t317 + t77;
t318 = t201 * t237;
t78 = t237 * t139 + t240 * t158;
t61 = -pkin(9) * t318 + t78;
t19 = -t245 * t61 + t248 * t56;
t408 = qJD(6) * t19 + t245 * t418 - t417 * t248;
t20 = t245 * t56 + t248 * t61;
t407 = -qJD(6) * t20 + t417 * t245 + t248 * t418;
t406 = t40 * mrSges(6,1);
t405 = t41 * mrSges(6,2);
t401 = t276 * Ifges(6,6);
t402 = t163 * Ifges(6,5);
t404 = t89 * Ifges(7,5) + Ifges(7,6) * t411 - t189 * Ifges(6,3) + t182 * Ifges(7,3) + t401 + t402;
t403 = mrSges(4,3) * t413;
t149 = t246 * t152;
t83 = t151 * t344 - t149;
t76 = -qJD(4) * pkin(4) + qJD(5) - t83;
t400 = t76 * t272;
t208 = t339 * t237;
t210 = t339 * t240;
t157 = -t208 * t245 + t210 * t248;
t321 = t189 * t240;
t135 = pkin(4) * t190 - qJ(5) * t189;
t52 = t240 * t135 - t237 * t83;
t35 = pkin(5) * t190 - pkin(9) * t321 + t52;
t322 = t189 * t237;
t53 = t237 * t135 + t240 * t83;
t42 = -pkin(9) * t322 + t53;
t398 = -qJD(5) * t200 - qJD(6) * t157 + t245 * t42 - t248 * t35;
t155 = -t208 * t248 - t210 * t245;
t397 = -qJD(5) * t263 + qJD(6) * t155 - t245 * t35 - t248 * t42;
t299 = qJDD(2) * t241;
t197 = -mrSges(4,1) * t299 + mrSges(4,2) * t300;
t72 = t144 * mrSges(5,1) + t143 * mrSges(5,2);
t396 = t197 + t72;
t58 = -t114 * mrSges(6,1) + t115 * mrSges(6,2);
t393 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t143 + t58;
t337 = mrSges(5,3) * t190;
t392 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t276 + mrSges(6,2) * t163 + t337;
t391 = Ifges(5,5) * qJD(4);
t390 = Ifges(5,6) * qJD(4);
t264 = t167 * t241 - (-t205 * t238 + t220) * t238;
t387 = t249 * t264;
t275 = -t241 * mrSges(4,1) + t238 * mrSges(4,2);
t382 = mrSges(5,1) * t189 - mrSges(5,2) * t190 - t275 * qJD(2);
t283 = qJD(2) * t309;
t215 = t249 * t283;
t303 = qJDD(1) * t239;
t181 = t247 * t303 + t215;
t165 = t181 + t414;
t302 = qJDD(1) * t242;
t218 = t241 * t302;
t133 = -t165 * t238 + t218;
t134 = t241 * t165 + t238 * t302;
t379 = -t133 * t238 + t134 * t241;
t109 = mrSges(6,2) * t189 + mrSges(6,3) * t276;
t110 = -mrSges(6,1) * t189 - mrSges(6,3) * t163;
t378 = t109 * t240 - t110 * t237;
t116 = t218 + (-pkin(8) * qJDD(2) - t165) * t238;
t117 = pkin(8) * t299 + t134;
t284 = qJD(4) * t344;
t293 = t246 * t116 + t344 * t117 + t151 * t284;
t25 = qJDD(4) * qJ(5) + (qJD(5) - t149) * qJD(4) + t293;
t214 = t247 * t283;
t180 = t249 * t303 - t214;
t262 = qJDD(3) - t180;
t153 = -qJDD(2) * t227 + t262;
t49 = pkin(4) * t144 - qJ(5) * t143 - qJD(5) * t190 + t153;
t12 = -t237 * t25 + t240 * t49;
t13 = t237 * t49 + t240 * t25;
t377 = -t12 * t237 + t13 * t240;
t376 = m(7) + m(6) + m(5);
t372 = -m(6) * t76 - t392;
t10 = pkin(9) * t114 + t13;
t5 = pkin(5) * t144 - pkin(9) * t115 + t12;
t1 = qJD(6) * t8 + t10 * t248 + t245 * t5;
t2 = -qJD(6) * t9 - t10 * t245 + t248 * t5;
t371 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t287 = m(4) * qJ(3) + mrSges(4,3);
t370 = mrSges(3,2) - t287 + t412;
t259 = m(4) * pkin(2) - t275;
t369 = mrSges(3,1) + t259 - t420;
t368 = Ifges(7,4) * t366 + Ifges(7,2) * t365 + Ifges(7,6) * t355;
t367 = Ifges(7,1) * t366 + Ifges(7,4) * t365 + Ifges(7,5) * t355;
t362 = Ifges(6,1) * t356 + Ifges(6,4) * t357 + Ifges(6,5) * t354;
t361 = -t411 / 0.2e1;
t360 = t411 / 0.2e1;
t359 = -t89 / 0.2e1;
t358 = t89 / 0.2e1;
t353 = -t182 / 0.2e1;
t352 = t182 / 0.2e1;
t351 = t189 / 0.2e1;
t350 = -t189 / 0.2e1;
t348 = t190 / 0.2e1;
t345 = t240 / 0.2e1;
t341 = g(3) * t239;
t338 = mrSges(5,3) * t189;
t336 = Ifges(5,4) * t190;
t335 = Ifges(6,4) * t237;
t334 = Ifges(6,4) * t240;
t330 = cos(pkin(10));
t329 = sin(pkin(10));
t316 = t239 * t247;
t307 = qJD(2) * t247;
t305 = qJD(4) * t246;
t298 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t138;
t48 = -mrSges(7,1) * t411 + mrSges(7,2) * t89;
t297 = t48 + t392;
t296 = m(4) + t376;
t295 = -t237 * (t163 * Ifges(6,4) + Ifges(6,2) * t276 - t189 * Ifges(6,6)) / 0.2e1;
t294 = (t163 * Ifges(6,1) + Ifges(6,4) * t276 - t189 * Ifges(6,5)) * t345;
t288 = t239 * t307;
t11 = -t34 * mrSges(7,1) + t33 * mrSges(7,2);
t282 = t239 * t330;
t281 = t239 * t329;
t280 = t330 * t247;
t279 = t330 * t249;
t278 = t329 * t247;
t277 = t329 * t249;
t270 = Ifges(6,1) * t240 - t335;
t269 = -Ifges(6,2) * t237 + t334;
t268 = Ifges(6,5) * t240 - Ifges(6,6) * t237;
t266 = t237 * t40 - t240 * t41;
t183 = -t238 * t316 + t241 * t242;
t184 = t238 * t242 + t241 * t316;
t121 = t246 * t183 + t184 * t344;
t97 = -t121 * t237 - t240 * t315;
t98 = t121 * t240 - t237 * t315;
t50 = -t245 * t98 + t248 * t97;
t51 = t245 * t97 + t248 * t98;
t261 = t183 * t344 - t246 * t184;
t28 = t116 * t344 - t246 * t117 - t151 * t305 - t152 * t284;
t186 = t242 * t280 + t277;
t145 = t186 * t230 + t232 * t282;
t188 = -t242 * t278 + t279;
t147 = t188 * t230 - t232 * t281;
t173 = t230 * t316 - t242 * t232;
t257 = -g(1) * t147 - g(2) * t145 - g(3) * t173;
t254 = t201 * t315;
t26 = -qJDD(4) * pkin(4) + qJDD(5) - t28;
t108 = qJD(3) * t201 + qJD(4) * t158;
t250 = qJD(2) ^ 2;
t203 = -qJD(2) * pkin(2) + t267;
t187 = t242 * t277 + t280;
t185 = -t242 * t279 + t278;
t179 = Ifges(5,4) * t189;
t174 = t230 * t242 + t232 * t316;
t171 = -qJD(4) * mrSges(5,2) + t338;
t170 = -qJDD(2) * pkin(2) + t262;
t148 = t188 * t232 + t230 * t281;
t146 = t186 * t232 - t230 * t282;
t131 = t190 * Ifges(5,1) + t179 + t391;
t130 = t189 * Ifges(5,2) + t336 + t390;
t127 = t263 * t201;
t126 = t200 * t201;
t123 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t144;
t111 = pkin(5) * t318 - t383;
t85 = pkin(5) * t320 + t108;
t82 = qJD(2) * t254 + qJD(4) * t121;
t81 = qJD(2) * t253 + qJD(4) * t261;
t70 = mrSges(7,1) * t182 - mrSges(7,3) * t89;
t69 = -mrSges(7,2) * t182 + mrSges(7,3) * t411;
t68 = mrSges(6,1) * t144 - mrSges(6,3) * t115;
t67 = -mrSges(6,2) * t144 + mrSges(6,3) * t114;
t66 = t237 * t288 + t240 * t81;
t65 = -t237 * t81 + t240 * t288;
t64 = t191 * t201 - t193 * t200;
t63 = -t192 * t201 - t193 * t263;
t62 = pkin(5) * t322 + t84;
t57 = -pkin(5) * t276 + t76;
t43 = t115 * Ifges(6,4) + t114 * Ifges(6,2) + t144 * Ifges(6,6);
t27 = -t152 * t305 + t293;
t22 = -mrSges(7,2) * t138 + mrSges(7,3) * t34;
t21 = mrSges(7,1) * t138 - mrSges(7,3) * t33;
t18 = -t114 * pkin(5) + t26;
t15 = -qJD(6) * t51 - t245 * t66 + t248 * t65;
t14 = qJD(6) * t50 + t245 * t65 + t248 * t66;
t3 = [t66 * t109 + t65 * t110 + t121 * t123 + t14 * t69 + t15 * t70 + t81 * t171 + t50 * t21 + t51 * t22 + t98 * t67 + t97 * t68 + (-t183 * t238 + t184 * t241) * qJDD(2) * mrSges(4,3) + t297 * t82 - (t11 + t393) * t261 + (-m(2) - m(3) - t296) * g(3) + m(6) * (t12 * t97 + t13 * t98 - t26 * t261 + t40 * t65 + t41 * t66 + t76 * t82) + m(7) * (t1 * t51 + t14 * t9 + t15 * t8 - t18 * t261 + t2 * t50 + t57 * t82) + m(5) * (t121 * t27 + t261 * t28 + t81 * t84 - t82 * t83) + m(4) * (t133 * t183 + t134 * t184) + ((mrSges(3,1) * qJDD(2) - t396) * t249 + (-qJDD(2) * mrSges(3,2) - qJD(2) * t382) * t247 + m(5) * (-t153 * t249 + t178 * t307) + m(4) * (qJD(2) * t387 - t170 * t249 + t203 * t307) + m(3) * (t180 * t249 + t181 * t247) + (-t247 * mrSges(3,1) + (-mrSges(3,2) + t403) * t249) * t250) * t239 + (m(3) * t242 ^ 2 + m(2)) * qJDD(1); (t153 * mrSges(5,2) - t28 * mrSges(5,3) + Ifges(5,1) * t143 - Ifges(5,4) * t144 + Ifges(5,5) * qJDD(4) + t26 * t272 + t268 * t354 + t269 * t357 + t270 * t356) * t201 + t392 * t108 + t394 * t109 + t395 * t110 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t360 + (-t12 * t317 - t13 * t318 - t319 * t40 - t320 * t41) * mrSges(6,3) + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t352 + (t294 + t295 + t400 + Ifges(5,1) * t348 + t268 * t350 + Ifges(5,4) * t351 + t276 * t269 / 0.2e1 + t163 * t270 / 0.2e1 + t131 / 0.2e1 + t391 / 0.2e1 + t178 * mrSges(5,2) - t83 * mrSges(5,3)) * t193 + (-pkin(2) * t170 + t264 * qJD(3) + t379 * qJ(3) - (t203 * t247 + t387) * t309) * m(4) - (Ifges(6,5) * t115 + Ifges(6,6) * t114 + Ifges(6,3) * t144 + t298) * t260 / 0.2e1 - (t153 * mrSges(5,1) + t12 * mrSges(6,1) - t13 * mrSges(6,2) - t27 * mrSges(5,3) - Ifges(5,4) * t143 + Ifges(6,5) * t356 + Ifges(7,5) * t366 + Ifges(5,2) * t144 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t357 + Ifges(7,6) * t365 + Ifges(6,3) * t354 + Ifges(7,3) * t355 + t371) * t260 + (-t376 * (-t185 * t227 + t186 * t340) + t370 * t186 + t369 * t185) * g(2) + (-t376 * (-t187 * t227 + t188 * t340) + t370 * t188 + t369 * t187) * g(1) + (t413 * (-t215 + t414) + t379) * mrSges(4,3) - t393 * t383 + (t108 * t76 + t12 * t77 + t13 * t78 - t26 * t383 + t394 * t41 + t395 * t40) * m(6) + (-t108 * t83 - t153 * t227 + t158 * t27 - t178 * t289 + t28 * t383 + t386 * t84) * m(5) + t382 * t289 + t386 * t171 + (t247 * t341 - t181 + t215) * mrSges(3,2) + (-Ifges(7,4) * t127 - Ifges(7,2) * t126) * t365 + (-Ifges(7,5) * t127 - Ifges(7,6) * t126) * t355 + (-t1 * t126 + t127 * t2 - t63 * t8 + t64 * t9) * mrSges(7,3) + (-Ifges(7,1) * t127 - Ifges(7,4) * t126) * t366 + t18 * (mrSges(7,1) * t126 - mrSges(7,2) * t127) + (m(5) * t83 - m(7) * t57 + t372 - t48) * qJD(1) * t254 + ((t420 * t249 + (-t340 * t376 + t412) * t247) * t239 - t376 * t227 * t315) * g(3) + (-t84 * mrSges(5,3) + t404 / 0.2e1 + t406 - t405 - Ifges(5,4) * t348 + Ifges(6,3) * t350 - Ifges(5,2) * t351 + Ifges(7,3) * t352 + Ifges(7,5) * t358 + Ifges(7,6) * t360 + t401 / 0.2e1 + t402 / 0.2e1 + t410 - t409 - t130 / 0.2e1 - t390 / 0.2e1 + t178 * mrSges(5,1)) * t194 + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t358 + Ifges(3,3) * qJDD(2) + (-t249 * t341 + t180 + t214) * mrSges(3,1) - (t247 * t287 + t249 * t259) * t341 - t127 * t367 - t126 * t368 - t43 * t318 / 0.2e1 + (Ifges(4,4) * t238 + Ifges(4,2) * t241) * t299 + (Ifges(4,1) * t238 + Ifges(4,4) * t241) * t300 + t407 * t70 + t408 * t69 + (t1 * t20 + t111 * t18 + t19 * t2 + t407 * t8 + t408 * t9 + t57 * t85) * m(7) + t19 * t21 + t20 * t22 + t317 * t362 + t63 * t363 + t64 * t364 + t57 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t77 * t68 + t78 * t67 + t85 * t48 + t111 * t11 + t158 * t123 + t170 * t275 - pkin(2) * t197 - t227 * t72; -t250 * t403 + t422 * t70 + t421 * t69 - t297 * t190 - (t171 + t378) * t189 + t240 * t68 - t263 * t21 + t200 * t22 + t237 * t67 + (t1 * t200 - t190 * t57 - t2 * t263 + t421 * t9 + t422 * t8) * m(7) + (t12 * t240 + t13 * t237 + t189 * t266 - t190 * t76) * m(6) + (-t189 * t84 + t190 * t83 + t153) * m(5) + (-qJD(2) * t264 + t170) * m(4) + (-g(1) * t187 - g(2) * t185 + g(3) * t315) * t296 + t396; t397 * t69 + (t1 * t157 + t155 * t2 - t18 * t226 + t397 * t9 + t398 * t8 - t57 * t62) * m(7) + t398 * t70 + (-Ifges(7,5) * t113 - Ifges(7,6) * t112 + Ifges(7,3) * t190) * t353 + (-Ifges(7,1) * t113 - Ifges(7,4) * t112 + Ifges(7,5) * t190) * t359 + (-Ifges(7,4) * t113 - Ifges(7,2) * t112 + Ifges(7,6) * t190) * t361 + (-Ifges(5,2) * t190 + t131 + t179) * t350 + t421 * t363 + t422 * t364 - t276 * (Ifges(6,6) * t190 + t189 * t269) / 0.2e1 - (Ifges(5,1) * t189 - t336 + t404) * t190 / 0.2e1 + (Ifges(6,3) * t190 + t189 * t268) * t351 - t163 * (Ifges(6,5) * t190 + t189 * t270) / 0.2e1 - qJD(4) * (Ifges(5,5) * t189 - Ifges(5,6) * t190) / 0.2e1 - t178 * (mrSges(5,1) * t190 + mrSges(5,2) * t189) + (-t1 * t263 - t2 * t200 - t421 * t8 + t422 * t9) * mrSges(7,3) + (-mrSges(7,1) * t422 + mrSges(7,2) * t421) * t57 + (t321 * t40 + t322 * t41 + t377) * mrSges(6,3) + (t338 - t171) * t83 + (t145 * t375 + t146 * t374) * g(2) + (t147 * t375 + t148 * t374) * g(1) + (t173 * t375 + t174 * t374) * g(3) + (Ifges(7,1) * t200 - Ifges(7,4) * t263) * t366 + t18 * (mrSges(7,1) * t263 + mrSges(7,2) * t200) + (Ifges(7,5) * t200 - Ifges(7,6) * t263) * t355 + (Ifges(7,4) * t200 - Ifges(7,2) * t263) * t365 - t263 * t368 + (t337 + t372) * t84 + (-t237 * t68 + t240 * t67) * qJ(5) + (-Ifges(7,5) * t191 - Ifges(7,6) * t192) * t352 + (-Ifges(7,1) * t191 - Ifges(7,4) * t192) * t358 + (-Ifges(7,4) * t191 - Ifges(7,2) * t192) * t360 + (-pkin(4) * t26 + qJ(5) * t377 - t266 * qJD(5) - t40 * t52 - t41 * t53) * m(6) + t378 * qJD(5) - t189 * t294 - t189 * t295 + Ifges(5,3) * qJDD(4) + t190 * t409 - t189 * t400 + t190 * t405 - t190 * t406 - t190 * t410 - t27 * mrSges(5,2) + t28 * mrSges(5,1) - pkin(4) * t58 - t62 * t48 + t43 * t345 + t130 * t348 + (Ifges(6,5) * t237 + Ifges(6,6) * t240) * t354 + (Ifges(6,1) * t237 + t334) * t356 + (Ifges(6,2) * t240 + t335) * t357 + t237 * t362 + t200 * t367 - t53 * t109 - t52 * t110 + Ifges(5,5) * t143 - Ifges(5,6) * t144 + t155 * t21 + t157 * t22 + t26 * t273 - t226 * t11; -t276 * t109 + t163 * t110 - t411 * t69 + t89 * t70 + t11 + t58 + (-t411 * t9 + t8 * t89 + t18 + t257) * m(7) + (t163 * t40 - t276 * t41 + t257 + t26) * m(6); -t57 * (mrSges(7,1) * t89 + mrSges(7,2) * t411) + (Ifges(7,1) * t411 - t343) * t359 + t38 * t358 + (Ifges(7,5) * t411 - Ifges(7,6) * t89) * t353 - t8 * t69 + t9 * t70 - g(1) * ((-t148 * t229 + t187 * t231) * mrSges(7,1) + (-t148 * t231 - t187 * t229) * mrSges(7,2)) - g(2) * ((-t146 * t229 + t185 * t231) * mrSges(7,1) + (-t146 * t231 - t185 * t229) * mrSges(7,2)) - g(3) * ((-t174 * t229 - t231 * t315) * mrSges(7,1) + (-t174 * t231 + t229 * t315) * mrSges(7,2)) + (t411 * t8 + t89 * t9) * mrSges(7,3) + t298 + (-Ifges(7,2) * t89 + t39 + t86) * t361 + t371;];
tau  = t3;

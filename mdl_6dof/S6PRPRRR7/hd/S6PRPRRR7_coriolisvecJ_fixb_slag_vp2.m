% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:07:07
% EndTime: 2018-11-23 15:07:23
% DurationCPUTime: 16.98s
% Computational Cost: add. (18839->657), mult. (60207->979), div. (0->0), fcn. (52256->16), ass. (0->316)
t237 = sin(pkin(14));
t240 = sin(pkin(6));
t241 = cos(pkin(14));
t252 = cos(qJ(2));
t243 = cos(pkin(7));
t248 = sin(qJ(2));
t330 = t243 * t248;
t259 = (-t237 * t252 - t241 * t330) * t240;
t196 = qJD(1) * t259;
t239 = sin(pkin(7));
t320 = qJD(3) * t239;
t429 = -t237 * t320 - t196;
t238 = sin(pkin(8));
t242 = cos(pkin(8));
t323 = qJD(1) * t240;
t309 = t248 * t323;
t302 = t239 * t309;
t401 = -t238 * t429 - t242 * t302;
t336 = t239 * t242;
t263 = (t238 * t243 + t241 * t336) * pkin(10);
t337 = t239 * t241;
t341 = t237 * t243;
t324 = pkin(2) * t341 + qJ(3) * t337;
t181 = t263 + t324;
t247 = sin(qJ(4));
t251 = cos(qJ(4));
t333 = t241 * t243;
t235 = pkin(2) * t333;
t342 = t237 * t239;
t381 = pkin(3) * t243;
t188 = t381 + t235 + (-pkin(10) * t242 - qJ(3)) * t342;
t380 = pkin(10) * t237;
t274 = -pkin(3) * t241 - t238 * t380;
t202 = (-pkin(2) + t274) * t239;
t277 = t188 * t242 + t202 * t238;
t112 = -t247 * t181 + t277 * t251;
t262 = t240 * (-t237 * t330 + t241 * t252);
t199 = qJD(1) * t262;
t332 = t242 * t247;
t260 = (-t237 * t332 + t241 * t251) * t239;
t282 = t238 * t302;
t406 = qJD(3) * t260 + qJD(4) * t112 - t199 * t251 - (t196 * t242 + t282) * t247;
t331 = t242 * t251;
t338 = t238 * t251;
t258 = t243 * t338 + t239 * (-t237 * t247 + t241 * t331);
t177 = t258 * qJD(4);
t339 = t238 * t247;
t185 = t243 * t339 + (t237 * t251 + t241 * t332) * t239;
t178 = t185 * qJD(4);
t428 = pkin(4) * t178 - pkin(11) * t177 + t401;
t175 = t258 * qJD(2);
t174 = qJD(5) - t175;
t201 = qJD(2) * t260;
t427 = -qJD(4) * t338 + t201;
t172 = t251 * t181;
t113 = t188 * t332 + t202 * t339 + t172;
t213 = -t238 * t337 + t242 * t243;
t100 = pkin(11) * t213 + t113;
t246 = sin(qJ(5));
t250 = cos(qJ(5));
t317 = qJD(5) * t250;
t318 = qJD(5) * t246;
t143 = -t188 * t238 + t242 * t202;
t95 = -pkin(4) * t258 - pkin(11) * t185 + t143;
t407 = -t100 * t318 + t428 * t246 + t250 * t406 + t95 * t317;
t261 = (t237 * t331 + t241 * t247) * t239;
t405 = qJD(3) * t261 + (t247 * t277 + t172) * qJD(4) + t196 * t331 - t199 * t247 + t251 * t282;
t321 = qJD(2) * t239;
t223 = qJ(3) * t321 + t309;
t308 = t252 * t323;
t227 = qJD(2) * pkin(2) + t308;
t244 = cos(pkin(6));
t322 = qJD(1) * t244;
t310 = t239 * t322;
t153 = t241 * t223 + t227 * t341 + t237 * t310;
t137 = qJD(2) * t263 + t153;
t152 = -t223 * t237 + t227 * t333 + t241 * t310;
t138 = (-t336 * t380 + t381) * qJD(2) + t152;
t314 = t243 * t322 + qJD(3);
t158 = (qJD(2) * t274 - t227) * t239 + t314;
t280 = t138 * t242 + t158 * t238;
t73 = t137 * t251 + t247 * t280;
t426 = -pkin(11) * qJD(6) * t250 - t73 + t174 * (pkin(5) * t246 - pkin(12) * t250);
t176 = t185 * qJD(2);
t230 = -pkin(5) * t250 - pkin(12) * t246 - pkin(4);
t131 = pkin(4) * t176 - pkin(11) * t175;
t72 = -t247 * t137 + t280 * t251;
t49 = t246 * t131 + t250 * t72;
t425 = pkin(11) * t318 + pkin(12) * t176 - qJD(6) * t230 + t49;
t424 = pkin(12) * t178 + t407;
t148 = t185 * t246 - t250 * t213;
t116 = -qJD(5) * t148 + t177 * t250;
t149 = t185 * t250 + t213 * t246;
t117 = qJD(5) * t149 + t177 * t246;
t423 = pkin(5) * t117 - pkin(12) * t116 + t405;
t220 = -t250 * t242 + t246 * t339;
t307 = t237 * t321;
t301 = t238 * t307;
t422 = qJD(5) * t220 + t246 * t301 + t250 * t427;
t198 = qJD(2) * t261;
t319 = qJD(4) * t247;
t421 = t238 * t319 - t198;
t204 = qJD(2) * t213 + qJD(4);
t67 = pkin(11) * t204 + t73;
t101 = -t138 * t238 + t242 * t158;
t74 = -pkin(4) * t175 - pkin(11) * t176 + t101;
t30 = -t246 * t67 + t250 * t74;
t28 = -pkin(5) * t174 - t30;
t245 = sin(qJ(6));
t249 = cos(qJ(6));
t31 = t246 * t74 + t250 * t67;
t29 = pkin(12) * t174 + t31;
t141 = -t176 * t246 + t204 * t250;
t142 = t176 * t250 + t204 * t246;
t66 = -pkin(4) * t204 - t72;
t36 = -pkin(5) * t141 - pkin(12) * t142 + t66;
t11 = -t245 * t29 + t249 * t36;
t12 = t245 * t36 + t249 * t29;
t287 = t11 * t249 + t12 * t245;
t289 = Ifges(7,5) * t249 - Ifges(7,6) * t245;
t366 = Ifges(7,4) * t249;
t291 = -Ifges(7,2) * t245 + t366;
t367 = Ifges(7,4) * t245;
t293 = Ifges(7,1) * t249 - t367;
t294 = mrSges(7,1) * t245 + mrSges(7,2) * t249;
t382 = t249 / 0.2e1;
t383 = -t245 / 0.2e1;
t140 = qJD(6) - t141;
t384 = t140 / 0.2e1;
t109 = t142 * t249 + t174 * t245;
t387 = t109 / 0.2e1;
t108 = -t142 * t245 + t174 * t249;
t389 = t108 / 0.2e1;
t368 = Ifges(7,4) * t109;
t46 = Ifges(7,2) * t108 + Ifges(7,6) * t140 + t368;
t105 = Ifges(7,4) * t108;
t47 = Ifges(7,1) * t109 + Ifges(7,5) * t140 + t105;
t420 = -t287 * mrSges(7,3) + t28 * t294 + t289 * t384 + t291 * t389 + t293 * t387 + t382 * t47 + t383 * t46;
t419 = Ifges(6,5) / 0.2e1;
t167 = qJD(2) * t177;
t106 = qJD(5) * t141 + t167 * t250;
t418 = t106 / 0.2e1;
t107 = qJD(5) * t142 + t167 * t246;
t417 = -t107 / 0.2e1;
t372 = t250 * t100 + t246 * t95;
t408 = -qJD(5) * t372 - t406 * t246 + t250 * t428;
t416 = t245 * t425 + t249 * t426;
t415 = t245 * t426 - t249 * t425;
t352 = t174 * Ifges(6,5);
t139 = Ifges(6,4) * t141;
t356 = t142 * Ifges(6,1);
t81 = t139 + t352 + t356;
t264 = -t81 / 0.2e1 - t352 / 0.2e1 - t66 * mrSges(6,2);
t414 = t30 * mrSges(6,3) + t264 - t420;
t413 = -t148 / 0.2e1;
t412 = t213 / 0.2e1;
t44 = -pkin(12) * t258 + t372;
t99 = -pkin(4) * t213 - t112;
t61 = pkin(5) * t148 - pkin(12) * t149 + t99;
t21 = -t245 * t44 + t249 * t61;
t411 = qJD(6) * t21 + t423 * t245 + t249 * t424;
t22 = t245 * t61 + t249 * t44;
t410 = -qJD(6) * t22 - t245 * t424 + t423 * t249;
t409 = -pkin(5) * t178 - t408;
t111 = mrSges(6,1) * t174 - mrSges(6,3) * t142;
t65 = -mrSges(7,1) * t108 + mrSges(7,2) * t109;
t346 = -t65 + t111;
t371 = mrSges(5,3) * t176;
t345 = -mrSges(5,1) * t204 - mrSges(6,1) * t141 + mrSges(6,2) * t142 + t371;
t221 = t242 * t246 + t250 * t339;
t271 = -t221 * t249 + t245 * t338;
t403 = qJD(6) * t271 + t245 * t422 + t249 * t421;
t192 = -t221 * t245 - t249 * t338;
t402 = qJD(6) * t192 + t245 * t421 - t249 * t422;
t400 = -qJD(5) * t221 + t246 * t427 - t250 * t301;
t215 = (t308 + t320) * qJD(2);
t334 = t240 * t248;
t306 = qJD(2) * t334;
t298 = qJD(1) * t306;
t275 = t243 * t298;
t179 = -t215 * t237 - t241 * t275;
t180 = t215 * t241 - t237 * t275;
t276 = t239 * t298;
t268 = t238 * t276;
t41 = t180 * t251 + (t179 * t242 + t268) * t247 + t72 * qJD(4);
t150 = -t179 * t238 + t242 * t276;
t168 = qJD(2) * t178;
t87 = pkin(4) * t168 - pkin(11) * t167 + t150;
t7 = t246 * t87 + t250 * t41 + t74 * t317 - t318 * t67;
t8 = -qJD(5) * t31 - t246 * t41 + t250 * t87;
t398 = -t246 * t8 + t250 * t7;
t42 = qJD(4) * t73 - t179 * t331 + t180 * t247 - t251 * t268;
t25 = pkin(5) * t107 - pkin(12) * t106 + t42;
t5 = pkin(12) * t168 + t7;
t1 = qJD(6) * t11 + t245 * t25 + t249 * t5;
t2 = -qJD(6) * t12 - t245 * t5 + t249 * t25;
t397 = t1 * t249 - t2 * t245;
t392 = Ifges(6,1) * t418 + Ifges(6,4) * t417 + t168 * t419;
t53 = qJD(6) * t108 + t106 * t249 + t168 * t245;
t54 = -qJD(6) * t109 - t106 * t245 + t168 * t249;
t15 = t53 * Ifges(7,1) + t54 * Ifges(7,4) + t107 * Ifges(7,5);
t396 = t15 / 0.2e1;
t395 = -t46 / 0.2e1;
t394 = t53 / 0.2e1;
t393 = t54 / 0.2e1;
t391 = t107 / 0.2e1;
t390 = -t108 / 0.2e1;
t388 = -t109 / 0.2e1;
t386 = -t139 / 0.2e1;
t385 = -t140 / 0.2e1;
t52 = Ifges(7,5) * t53;
t51 = Ifges(7,6) * t54;
t375 = t31 * mrSges(6,3);
t374 = t72 * mrSges(5,3);
t20 = -mrSges(7,1) * t54 + mrSges(7,2) * t53;
t84 = mrSges(6,1) * t168 - mrSges(6,3) * t106;
t373 = t20 - t84;
t370 = Ifges(5,4) * t176;
t369 = Ifges(6,4) * t142;
t364 = t106 * Ifges(6,4);
t362 = t108 * Ifges(7,6);
t361 = t109 * Ifges(7,5);
t329 = t243 * t252;
t335 = t239 * t244;
t182 = t241 * t335 + (-t237 * t248 + t241 * t329) * t240;
t214 = -t239 * t240 * t252 + t243 * t244;
t183 = t241 * t334 + (t240 * t329 + t335) * t237;
t344 = t183 * t247;
t118 = -t182 * t331 - t214 * t338 + t344;
t360 = t118 * t42;
t359 = t140 * Ifges(7,3);
t358 = t141 * Ifges(6,2);
t357 = t141 * Ifges(6,6);
t355 = t142 * Ifges(6,5);
t353 = t168 * Ifges(6,6);
t351 = t174 * Ifges(6,6);
t350 = t174 * Ifges(6,3);
t349 = t204 * Ifges(5,5);
t348 = t204 * Ifges(5,6);
t347 = t251 * t42;
t272 = -mrSges(4,2) * t243 + mrSges(4,3) * t337;
t219 = t272 * qJD(2);
t343 = t219 * t241;
t328 = t245 * t250;
t327 = t249 * t250;
t325 = Ifges(5,5) * t167 - Ifges(5,6) * t168;
t13 = Ifges(7,3) * t107 + t51 + t52;
t311 = Ifges(6,5) * t106 - Ifges(6,6) * t107 + Ifges(6,3) * t168;
t300 = t239 * t306;
t297 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t295 = mrSges(7,1) * t249 - mrSges(7,2) * t245;
t292 = Ifges(7,1) * t245 + t366;
t290 = Ifges(7,2) * t249 + t367;
t288 = Ifges(7,5) * t245 + Ifges(7,6) * t249;
t286 = t246 * t31 + t250 * t30;
t59 = -t100 * t246 + t250 * t95;
t278 = t182 * t242 + t214 * t238;
t119 = t183 * t251 + t247 * t278;
t147 = -t182 * t238 + t214 * t242;
t83 = t119 * t250 + t147 * t246;
t57 = t118 * t249 - t245 * t83;
t58 = t118 * t245 + t249 * t83;
t48 = t131 * t250 - t246 * t72;
t281 = t238 * t300;
t82 = t119 * t246 - t250 * t147;
t121 = t149 * t249 - t245 * t258;
t120 = -t149 * t245 - t249 * t258;
t279 = -t152 * t237 + t153 * t241;
t273 = mrSges(4,1) * t243 - mrSges(4,3) * t342;
t267 = qJD(2) * (-mrSges(4,1) * t241 + mrSges(4,2) * t237);
t45 = t359 + t361 + t362;
t80 = t351 + t358 + t369;
t257 = t12 * mrSges(7,2) - t45 / 0.2e1 + t80 / 0.2e1 + t369 / 0.2e1 - t362 / 0.2e1 - t361 / 0.2e1 - t11 * mrSges(7,1) - t359 / 0.2e1 + t351 / 0.2e1 - t66 * mrSges(6,1);
t256 = t358 / 0.2e1 + t257;
t253 = qJD(2) ^ 2;
t218 = t273 * qJD(2);
t208 = t239 * t267;
t206 = pkin(11) * t327 + t230 * t245;
t205 = -pkin(11) * t328 + t230 * t249;
t200 = qJD(2) * t262;
t197 = qJD(2) * t259;
t195 = -t227 * t239 + t314;
t173 = Ifges(5,4) * t175;
t160 = -t197 * t238 + t242 * t300;
t144 = -mrSges(5,2) * t204 + mrSges(5,3) * t175;
t130 = -mrSges(5,1) * t175 + mrSges(5,2) * t176;
t126 = t175 * t327 + t176 * t245;
t125 = -t175 * t328 + t176 * t249;
t124 = mrSges(5,1) * t168 + mrSges(5,2) * t167;
t115 = t176 * Ifges(5,1) + t173 + t349;
t114 = t175 * Ifges(5,2) + t348 + t370;
t110 = -mrSges(6,2) * t174 + mrSges(6,3) * t141;
t93 = pkin(5) * t142 - pkin(12) * t141;
t85 = -mrSges(6,2) * t168 - mrSges(6,3) * t107;
t79 = t350 + t355 + t357;
t78 = qJD(4) * t119 - t197 * t331 + t200 * t247 - t251 * t281;
t77 = t200 * t251 + (t197 * t242 + t281) * t247 + (t251 * t278 - t344) * qJD(4);
t76 = mrSges(7,1) * t140 - mrSges(7,3) * t109;
t75 = -mrSges(7,2) * t140 + mrSges(7,3) * t108;
t64 = mrSges(6,1) * t107 + mrSges(6,2) * t106;
t63 = -qJD(6) * t121 - t116 * t245 + t178 * t249;
t62 = qJD(6) * t120 + t116 * t249 + t178 * t245;
t55 = -t107 * Ifges(6,2) + t353 + t364;
t43 = pkin(5) * t258 - t59;
t38 = -pkin(5) * t176 - t48;
t35 = -mrSges(7,2) * t107 + mrSges(7,3) * t54;
t34 = mrSges(7,1) * t107 - mrSges(7,3) * t53;
t33 = qJD(5) * t83 - t250 * t160 + t246 * t77;
t32 = -qJD(5) * t82 + t160 * t246 + t250 * t77;
t19 = t245 * t93 + t249 * t30;
t18 = -t245 * t30 + t249 * t93;
t14 = t53 * Ifges(7,4) + t54 * Ifges(7,2) + t107 * Ifges(7,6);
t10 = -qJD(6) * t58 - t245 * t32 + t249 * t78;
t9 = qJD(6) * t57 + t245 * t78 + t249 * t32;
t6 = -pkin(5) * t168 - t8;
t3 = [t10 * t76 + t32 * t110 + t118 * t64 + t147 * t124 + t160 * t130 + t77 * t144 + t197 * t218 + t200 * t219 + t57 * t34 + t58 * t35 + t9 * t75 + t83 * t85 + t373 * t82 + t345 * t78 - t346 * t33 + (t118 * t167 - t119 * t168) * mrSges(5,3) + (-t253 * t252 * mrSges(3,2) + (-mrSges(3,1) * t253 + t208 * t321) * t248) * t240 + m(6) * (-t30 * t33 + t31 * t32 + t66 * t78 + t7 * t83 - t8 * t82 + t360) + m(7) * (t1 * t58 + t10 * t11 + t12 * t9 + t2 * t57 + t28 * t33 + t6 * t82) + m(5) * (t101 * t160 + t119 * t41 + t147 * t150 - t72 * t78 + t73 * t77 + t360) + m(4) * (t152 * t197 + t153 * t200 + t179 * t182 + t180 * t183 + (qJD(1) * t214 + t195) * t300); -t177 * t374 + t405 * t345 + t401 * t130 + (-mrSges(5,1) * t213 + mrSges(6,1) * t148 + mrSges(6,2) * t149 + mrSges(5,3) * t185) * t42 + t325 * t412 + t55 * t413 + t99 * t64 + t59 * t84 + (t180 * t324 + t179 * t235 + (-t179 * qJ(3) * t237 - pkin(2) * t276 + qJD(3) * t279) * t239 - t152 * t196 - t153 * t199 - t195 * t302) * m(4) + t180 * t272 + t179 * t273 - t258 * t311 / 0.2e1 + (-t113 * mrSges(5,3) + t149 * t419 + Ifges(6,6) * t413 - Ifges(5,4) * t185 - Ifges(5,6) * t213 / 0.2e1 - (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t258) * t168 + (-t112 * mrSges(5,3) + Ifges(5,1) * t185 + Ifges(5,4) * t258 + Ifges(5,5) * t412) * t167 + (Ifges(6,4) * t149 - Ifges(6,2) * t148 - Ifges(6,6) * t258) * t417 + (Ifges(6,1) * t149 - Ifges(6,4) * t148 - Ifges(6,5) * t258) * t418 + t7 * (mrSges(6,2) * t258 - mrSges(6,3) * t148) + t8 * (-mrSges(6,1) * t258 - mrSges(6,3) * t149) + t150 * (-mrSges(5,1) * t258 + mrSges(5,2) * t185) + t41 * (-mrSges(5,2) * t213 + mrSges(5,3) * t258) - t208 * t302 + t372 * t85 + t429 * t218 + t43 * t20 + t21 * t34 + t22 * t35 + (Ifges(7,4) * t62 + Ifges(7,2) * t63 + Ifges(7,6) * t117) * t389 + (Ifges(7,5) * t121 + Ifges(7,6) * t120 + Ifges(7,3) * t148) * t391 + t149 * t392 + (Ifges(7,4) * t121 + Ifges(7,2) * t120 + Ifges(7,6) * t148) * t393 + (Ifges(7,1) * t121 + Ifges(7,4) * t120 + Ifges(7,5) * t148) * t394 + t121 * t396 + (Ifges(7,5) * t62 + Ifges(7,6) * t63 + Ifges(7,3) * t117) * t384 + (Ifges(7,1) * t62 + Ifges(7,4) * t63 + Ifges(7,5) * t117) * t387 + t320 * t343 - t73 * t178 * mrSges(5,3) + t116 * t81 / 0.2e1 + t117 * t45 / 0.2e1 + t11 * (mrSges(7,1) * t117 - mrSges(7,3) * t62) + t12 * (-mrSges(7,2) * t117 + mrSges(7,3) * t63) - t117 * t80 / 0.2e1 + t66 * (mrSges(6,1) * t117 + mrSges(6,2) * t116) + t120 * t14 / 0.2e1 + t6 * (-mrSges(7,1) * t120 + mrSges(7,2) * t121) + t406 * t144 + (t101 * t401 - t112 * t42 + t113 * t41 + t143 * t150 - t405 * t72 + t406 * t73) * m(5) + t143 * t124 + t148 * t13 / 0.2e1 + t1 * (-mrSges(7,2) * t148 + mrSges(7,3) * t120) + t2 * (mrSges(7,1) * t148 - mrSges(7,3) * t121) + t407 * t110 + t408 * t111 + (t30 * t408 + t31 * t407 + t372 * t7 + t405 * t66 + t42 * t99 + t59 * t8) * m(6) + t409 * t65 + t410 * t76 + t411 * t75 + (t1 * t22 + t11 * t410 + t12 * t411 + t2 * t21 + t28 * t409 + t43 * t6) * m(7) + t177 * t115 / 0.2e1 + t62 * t47 / 0.2e1 + t63 * t46 / 0.2e1 + t28 * (-mrSges(7,1) * t63 + mrSges(7,2) * t62) + t178 * t79 / 0.2e1 + t30 * (mrSges(6,1) * t178 - mrSges(6,3) * t116) + t31 * (-mrSges(6,2) * t178 - mrSges(6,3) * t117) + t141 * (Ifges(6,4) * t116 - Ifges(6,2) * t117 + Ifges(6,6) * t178) / 0.2e1 + t142 * (Ifges(6,1) * t116 - Ifges(6,4) * t117 + Ifges(6,5) * t178) / 0.2e1 + t174 * (Ifges(6,5) * t116 - Ifges(6,6) * t117 + Ifges(6,3) * t178) / 0.2e1 - t178 * t114 / 0.2e1 + t101 * (mrSges(5,1) * t178 + mrSges(5,2) * t177) + t175 * (Ifges(5,4) * t177 - Ifges(5,2) * t178) / 0.2e1 + t176 * (Ifges(5,1) * t177 - Ifges(5,4) * t178) / 0.2e1 + t204 * (Ifges(5,5) * t177 - Ifges(5,6) * t178) / 0.2e1 - t199 * t219 + t239 ^ 2 * t267 * t309; t242 * t124 - t201 * t144 + t192 * t34 - t271 * t35 + t221 * t85 + t403 * t76 + t402 * t75 + t373 * t220 - t345 * t198 - t422 * t110 + (-t130 * t307 - t251 * t64 + (-t167 * t251 - t168 * t247) * mrSges(5,3) + (t144 * t251 + t247 * t345) * qJD(4)) * t238 + t346 * t400 + (t237 * t218 - t343 + (-t279 + t309) * m(4)) * t321 + (-t1 * t271 + t11 * t403 + t12 * t402 + t192 * t2 + t220 * t6 - t28 * t400) * m(7) + (-t198 * t66 - t220 * t8 + t221 * t7 + (t319 * t66 - t347) * t238 - t422 * t31 + t400 * t30) * m(6) + (-t101 * t301 + t198 * t72 - t201 * t73 + t150 * t242 + (t247 * t41 - t347 + (-t247 * t72 + t251 * t73) * qJD(4)) * t238) * m(5); (-t345 + t371) * t73 + t415 * t75 + (t1 * t206 + t11 * t416 + t12 * t415 + t2 * t205 - t28 * t38) * m(7) + t416 * t76 + t325 + (t175 * t286 + t398) * mrSges(6,3) + (t11 * t126 - t12 * t125) * mrSges(7,3) + (t85 * pkin(11) + t55 / 0.2e1 - t13 / 0.2e1 - t42 * mrSges(6,1) + t364 / 0.2e1 + t353 / 0.2e1 - t52 / 0.2e1 - t51 / 0.2e1 + (-Ifges(6,2) / 0.2e1 - Ifges(7,3) / 0.2e1) * t107 + t297 + (t386 - t356 / 0.2e1 + t264) * t175 + (t139 / 0.2e1 + t356 / 0.2e1 + (m(7) * t28 - t346) * pkin(11) - t414) * qJD(5)) * t250 + (t374 - t115 / 0.2e1 - t101 * mrSges(5,2) - t173 / 0.2e1 - t349 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t176) * t175 + (-t101 * mrSges(5,1) + t31 * mrSges(6,2) - t30 * mrSges(6,1) - t357 / 0.2e1 - t355 / 0.2e1 - t350 / 0.2e1 + t348 / 0.2e1 - t79 / 0.2e1 + t114 / 0.2e1 + t370 / 0.2e1) * t176 + ((t28 * t295 + t290 * t390 + t292 * t388 + t288 * t385 + t47 * t383 + t249 * t395 + (t11 * t245 - t12 * t249) * mrSges(7,3)) * qJD(6) + t256 * t175 + t42 * mrSges(6,2) + t14 * t383 + t15 * t382 + t289 * t391 + t291 * t393 + t293 * t394 + t6 * t294 + 0.2e1 * t392 + (-t256 - t375) * qJD(5) + (-t1 * t245 - t2 * t249) * mrSges(7,3) + (m(7) * t6 - qJD(5) * t110 + t373) * pkin(11)) * t246 + (-pkin(4) * t42 - t30 * t48 - t31 * t49 - t66 * t73 + (-qJD(5) * t286 + t398) * pkin(11)) * m(6) - t41 * mrSges(5,2) - t42 * mrSges(5,1) + (Ifges(7,4) * t126 + Ifges(7,2) * t125) * t390 + t125 * t395 + (Ifges(7,5) * t126 + Ifges(7,6) * t125) * t385 + (Ifges(7,1) * t126 + Ifges(7,4) * t125) * t388 - t49 * t110 - t48 * t111 - t126 * t47 / 0.2e1 - t28 * (-mrSges(7,1) * t125 + mrSges(7,2) * t126) - t72 * t144 - pkin(4) * t64 - t38 * t65 + t205 * t34 + t206 * t35; -t6 * t295 - t19 * t75 - t18 * t76 + t397 * mrSges(7,3) + (t257 + t375) * t142 + t346 * t31 + (t386 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t142 + t414) * t141 + t420 * qJD(6) + t288 * t391 + t290 * t393 + t292 * t394 - pkin(5) * t20 - t7 * mrSges(6,2) + t8 * mrSges(6,1) - t30 * t110 + t311 + t245 * t396 + t14 * t382 + (-pkin(5) * t6 - t11 * t18 - t12 * t19 - t28 * t31) * m(7) + (m(7) * t397 - t245 * t34 + t249 * t35 + (-m(7) * t287 - t245 * t75 - t249 * t76) * qJD(6)) * pkin(12); -t28 * (mrSges(7,1) * t109 + mrSges(7,2) * t108) + t12 * t76 + (Ifges(7,1) * t108 - t368) * t388 + t46 * t387 + (Ifges(7,5) * t108 - Ifges(7,6) * t109) * t385 - t11 * t75 + (t108 * t11 + t109 * t12) * mrSges(7,3) - t297 + t13 + (-Ifges(7,2) * t109 + t105 + t47) * t390;];
tauc  = t3(:);

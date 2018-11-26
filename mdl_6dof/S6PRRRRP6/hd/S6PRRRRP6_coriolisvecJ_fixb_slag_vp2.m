% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:31:10
% EndTime: 2018-11-23 15:31:28
% DurationCPUTime: 18.59s
% Computational Cost: add. (10842->720), mult. (29436->973), div. (0->0), fcn. (23236->12), ass. (0->321)
t244 = cos(pkin(7));
t252 = cos(qJ(3));
t253 = cos(qJ(2));
t346 = t252 * t253;
t248 = sin(qJ(3));
t249 = sin(qJ(2));
t350 = t248 * t249;
t271 = -t244 * t350 + t346;
t243 = sin(pkin(6));
t341 = qJD(1) * t243;
t174 = t271 * t341;
t242 = sin(pkin(7));
t356 = t242 * t248;
t236 = pkin(9) * t356;
t352 = t244 * t252;
t213 = pkin(2) * t352 - t236;
t205 = t213 * qJD(3);
t420 = t205 - t174;
t432 = Ifges(6,1) + Ifges(7,1);
t431 = Ifges(7,4) + Ifges(6,5);
t353 = t244 * t248;
t355 = t242 * t252;
t214 = pkin(2) * t353 + pkin(9) * t355;
t197 = pkin(10) * t244 + t214;
t305 = -pkin(3) * t252 - pkin(10) * t248;
t198 = (-pkin(2) + t305) * t242;
t304 = pkin(3) * t248 - pkin(10) * t252;
t268 = t304 * qJD(3);
t204 = t242 * t268;
t247 = sin(qJ(4));
t251 = cos(qJ(4));
t321 = t249 * t341;
t313 = t242 * t321;
t334 = qJD(4) * t251;
t335 = qJD(4) * t247;
t424 = -t197 * t335 + t198 * t334 + t420 * t251 + (t204 - t313) * t247;
t348 = t249 * t252;
t349 = t248 * t253;
t273 = t244 * t348 + t349;
t419 = t214 * qJD(3) - t273 * t341;
t337 = qJD(2) * t252;
t319 = t242 * t337;
t453 = t247 * t319 - t335;
t430 = Ifges(7,5) - Ifges(6,4);
t447 = Ifges(6,6) - Ifges(7,6);
t303 = pkin(4) * t247 - pkin(11) * t251;
t221 = t303 * qJD(4);
t231 = -pkin(4) * t251 - pkin(11) * t247 - pkin(3);
t246 = sin(qJ(5));
t250 = cos(qJ(5));
t331 = qJD(5) * t251;
t332 = qJD(5) * t250;
t135 = t246 * t221 + t231 * t332 + (-t246 * t331 - t250 * t335) * pkin(10);
t224 = qJD(2) * pkin(2) + t253 * t341;
t339 = qJD(2) * t242;
t216 = pkin(9) * t339 + t321;
t358 = t216 * t252;
t274 = t224 * t353 + t358;
t245 = cos(pkin(6));
t340 = qJD(1) * t245;
t103 = (t248 * t340 + t303 * t337) * t242 + t274;
t207 = t248 * t216;
t266 = t224 * t244 + t242 * t340;
t138 = t252 * t266 - t207;
t203 = t304 * t339;
t102 = t251 * t138 + t247 * t203;
t320 = t248 * t339;
t91 = pkin(11) * t320 + t102;
t38 = t246 * t103 + t250 * t91;
t452 = t135 - t38;
t333 = qJD(5) * t246;
t136 = pkin(10) * (t246 * t335 - t250 * t331) + t221 * t250 - t231 * t333;
t37 = t103 * t250 - t246 * t91;
t451 = t136 - t37;
t316 = qJD(3) * t356;
t450 = -pkin(11) * t316 - t424;
t210 = -t251 * t244 + t247 * t356;
t336 = qJD(3) * t252;
t317 = t251 * t336;
t164 = -qJD(4) * t210 + t242 * t317;
t211 = t244 * t247 + t251 * t356;
t318 = t247 * t336;
t165 = qJD(4) * t211 + t242 * t318;
t449 = pkin(4) * t165 - pkin(11) * t164 + t419;
t338 = qJD(2) * t244;
t235 = qJD(3) + t338;
t187 = t235 * t247 + t251 * t320;
t448 = -Ifges(5,4) * t187 / 0.2e1;
t282 = -qJD(4) + t319;
t149 = t246 * t187 + t250 * t282;
t148 = Ifges(6,4) * t149;
t150 = t250 * t187 - t246 * t282;
t186 = t235 * t251 - t247 * t320;
t182 = qJD(5) - t186;
t371 = Ifges(7,5) * t149;
t425 = t150 * t432 + t431 * t182 - t148 + t371;
t446 = pkin(5) * t453 - t451;
t351 = t246 * t251;
t178 = -t250 * t320 + t319 * t351;
t347 = t250 * t252;
t179 = (t246 * t248 + t251 * t347) * t339;
t283 = pkin(5) * t246 - qJ(6) * t250;
t275 = pkin(10) + t283;
t284 = pkin(5) * t250 + qJ(6) * t246;
t101 = -t247 * t138 + t203 * t251;
t90 = -pkin(4) * t320 - t101;
t445 = -pkin(5) * t178 + qJ(6) * t179 + (qJD(5) * t284 - qJD(6) * t250) * t247 + t275 * t334 - t90;
t444 = -qJ(6) * t453 - qJD(6) * t251 + t452;
t368 = Ifges(5,2) * t186;
t120 = -pkin(3) * t235 - t138;
t139 = t248 * t266 + t358;
t121 = pkin(10) * t235 + t139;
t234 = t244 * t340;
t146 = t234 + (qJD(2) * t305 - t224) * t242;
t65 = t251 * t121 + t247 * t146;
t54 = -pkin(11) * t282 + t65;
t67 = -pkin(4) * t186 - pkin(11) * t187 + t120;
t24 = -t246 * t54 + t250 * t67;
t422 = qJD(6) - t24;
t13 = -pkin(5) * t182 + t422;
t25 = t246 * t67 + t250 * t54;
t14 = qJ(6) * t182 + t25;
t325 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t326 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t327 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t269 = Ifges(5,6) * t282;
t441 = t269 / 0.2e1 - t368 / 0.2e1 + t448;
t59 = Ifges(6,5) * t150 - Ifges(6,6) * t149 + Ifges(6,3) * t182;
t60 = Ifges(7,4) * t150 + Ifges(7,2) * t182 + Ifges(7,6) * t149;
t410 = t325 * t149 - t327 * t150 + t326 * t182 - t13 * mrSges(7,1) - t25 * mrSges(6,2) - t65 * mrSges(5,3) + t441 + t59 / 0.2e1 + t60 / 0.2e1 + t448 + t120 * mrSges(5,1) + t14 * mrSges(7,3) + t24 * mrSges(6,1);
t443 = t368 / 0.2e1 - t410;
t64 = -t247 * t121 + t251 * t146;
t53 = pkin(4) * t282 - t64;
t28 = t149 * pkin(5) - t150 * qJ(6) + t53;
t147 = Ifges(7,5) * t150;
t58 = Ifges(7,6) * t182 + Ifges(7,3) * t149 + t147;
t375 = Ifges(6,4) * t150;
t61 = -Ifges(6,2) * t149 + Ifges(6,6) * t182 + t375;
t442 = t28 * mrSges(7,1) + t53 * mrSges(6,1) + t58 / 0.2e1 - t61 / 0.2e1 - t14 * mrSges(7,2) - t25 * mrSges(6,3);
t391 = t235 / 0.2e1;
t196 = t236 + (-pkin(2) * t252 - pkin(3)) * t244;
t124 = pkin(4) * t210 - pkin(11) * t211 + t196;
t141 = t251 * t197 + t247 * t198;
t126 = -pkin(11) * t355 + t141;
t434 = t124 * t332 - t126 * t333 + t246 * t449 - t250 * t450;
t421 = t246 * t124 + t250 * t126;
t433 = -qJD(5) * t421 + t246 * t450 + t250 * t449;
t156 = t235 * t335 + (t248 * t334 + t318) * t339;
t155 = t235 * t334 + (-t248 * t335 + t317) * t339;
t315 = qJD(3) * t339;
t308 = t248 * t315;
t77 = -qJD(5) * t149 + t250 * t155 + t246 * t308;
t78 = qJD(5) * t150 + t246 * t155 - t250 * t308;
t426 = t431 * t156 + t430 * t78 + t432 * t77;
t440 = -t246 * t447 + t250 * t431;
t370 = Ifges(7,5) * t246;
t374 = Ifges(6,4) * t246;
t439 = t250 * t432 + t370 - t374;
t438 = t320 / 0.2e1;
t437 = -t235 * Ifges(4,6) / 0.2e1;
t436 = qJ(6) * t165 + qJD(6) * t210 + t434;
t435 = -pkin(5) * t165 - t433;
t429 = Ifges(7,2) + Ifges(6,3);
t19 = Ifges(6,5) * t77 - Ifges(6,6) * t78 + Ifges(6,3) * t156;
t20 = Ifges(7,4) * t77 + Ifges(7,2) * t156 + Ifges(7,6) * t78;
t427 = t19 + t20;
t423 = -qJD(6) * t246 + t182 * t283 - t65;
t343 = mrSges(4,1) * t235 + mrSges(5,1) * t186 - mrSges(5,2) * t187 - mrSges(4,3) * t320;
t181 = Ifges(5,4) * t186;
t270 = Ifges(5,5) * t282;
t377 = Ifges(5,1) * t187;
t113 = t181 - t270 + t377;
t418 = t64 * mrSges(5,3) - t113 / 0.2e1 - t120 * mrSges(5,2) - t181 / 0.2e1;
t417 = t246 * t431 + t250 * t447;
t369 = Ifges(7,5) * t250;
t373 = Ifges(6,4) * t250;
t416 = t246 * t432 - t369 + t373;
t415 = t244 * t346 - t350;
t161 = (t268 + t321) * t339;
t262 = t271 * qJD(2);
t323 = t245 * t355;
t309 = qJD(3) * t323;
t93 = (t224 * t352 - t207) * qJD(3) + (t243 * t262 + t309) * qJD(1);
t26 = -t121 * t335 + t146 * t334 + t247 * t161 + t251 * t93;
t16 = pkin(11) * t308 + t26;
t263 = t273 * qJD(2);
t310 = t245 * t316;
t94 = t274 * qJD(3) + (t243 * t263 + t310) * qJD(1);
t41 = pkin(4) * t156 - pkin(11) * t155 + t94;
t3 = t250 * t16 + t246 * t41 + t67 * t332 - t333 * t54;
t4 = -qJD(5) * t25 - t16 * t246 + t250 * t41;
t414 = -t246 * t4 + t250 * t3;
t1 = qJ(6) * t156 + qJD(6) * t182 + t3;
t2 = -pkin(5) * t156 - t4;
t413 = t1 * t250 + t2 * t246;
t27 = -t121 * t334 - t146 * t335 + t161 * t251 - t247 * t93;
t412 = -t27 * mrSges(5,1) + t26 * mrSges(5,2) - Ifges(5,5) * t155 + Ifges(5,6) * t156;
t411 = Ifges(4,4) * t319 / 0.2e1 + Ifges(4,5) * t391;
t279 = t24 * t250 + t246 * t25;
t281 = t13 * t250 - t14 * t246;
t286 = Ifges(7,3) * t246 + t369;
t293 = -Ifges(6,2) * t246 + t373;
t298 = mrSges(7,1) * t246 - mrSges(7,3) * t250;
t300 = mrSges(6,1) * t246 + mrSges(6,2) * t250;
t387 = t250 / 0.2e1;
t389 = t246 / 0.2e1;
t390 = -t246 / 0.2e1;
t395 = t182 / 0.2e1;
t402 = t150 / 0.2e1;
t404 = t149 / 0.2e1;
t405 = -t149 / 0.2e1;
t255 = t281 * mrSges(7,2) - t279 * mrSges(6,3) + t28 * t298 + t286 * t404 + t293 * t405 + t300 * t53 + t425 * t387 + t389 * t58 + t390 * t61 + t395 * t440 + t402 * t439;
t409 = t77 / 0.2e1;
t408 = -t78 / 0.2e1;
t407 = t78 / 0.2e1;
t403 = -t150 / 0.2e1;
t401 = t156 / 0.2e1;
t396 = -t182 / 0.2e1;
t394 = -t210 / 0.2e1;
t392 = t211 / 0.2e1;
t388 = -t250 / 0.2e1;
t382 = qJD(4) / 0.2e1;
t46 = mrSges(6,1) * t156 - mrSges(6,3) * t77;
t47 = -t156 * mrSges(7,1) + t77 * mrSges(7,2);
t381 = t47 - t46;
t48 = -mrSges(6,2) * t156 - mrSges(6,3) * t78;
t49 = -mrSges(7,2) * t78 + mrSges(7,3) * t156;
t380 = t48 + t49;
t379 = mrSges(6,3) * t149;
t378 = mrSges(6,3) * t150;
t366 = t155 * Ifges(5,1);
t365 = t155 * Ifges(5,4);
t364 = t156 * Ifges(5,4);
t162 = -t243 * t415 - t323;
t363 = t162 * t94;
t131 = mrSges(5,1) * t308 - mrSges(5,3) * t155;
t30 = mrSges(6,1) * t78 + mrSges(6,2) * t77;
t360 = t30 - t131;
t137 = pkin(4) * t187 - pkin(11) * t186;
t36 = t246 * t137 + t250 * t64;
t159 = -mrSges(5,1) * t282 - t187 * mrSges(5,3);
t87 = mrSges(6,1) * t149 + mrSges(6,2) * t150;
t359 = t87 - t159;
t357 = t231 * t250;
t354 = t243 * t249;
t104 = -mrSges(7,2) * t149 + mrSges(7,3) * t182;
t105 = -mrSges(6,2) * t182 - t379;
t345 = t104 + t105;
t106 = mrSges(6,1) * t182 - t378;
t107 = -mrSges(7,1) * t182 + mrSges(7,2) * t150;
t344 = -t106 + t107;
t194 = t250 * t251 * pkin(10) + t246 * t231;
t86 = mrSges(7,1) * t149 - mrSges(7,3) * t150;
t328 = -t86 - t359;
t324 = t246 * t355;
t140 = -t247 * t197 + t198 * t251;
t180 = -t224 * t242 + t234;
t314 = t180 * mrSges(4,2) + Ifges(4,1) * t438 + t411;
t312 = t339 * t354;
t307 = -t94 * mrSges(4,1) - t93 * mrSges(4,2);
t125 = pkin(4) * t355 - t140;
t302 = -mrSges(4,1) * t252 + mrSges(4,2) * t248;
t301 = mrSges(6,1) * t250 - mrSges(6,2) * t246;
t299 = mrSges(7,1) * t250 + mrSges(7,3) * t246;
t292 = Ifges(6,2) * t250 + t374;
t289 = Ifges(5,5) * t164 - Ifges(5,6) * t165;
t285 = -Ifges(7,3) * t250 + t370;
t35 = t137 * t250 - t246 * t64;
t272 = t244 * t349 + t348;
t163 = t243 * t272 + t245 * t356;
t209 = -t242 * t243 * t253 + t244 * t245;
t123 = t163 * t251 + t209 * t247;
t76 = t123 * t250 + t162 * t246;
t75 = t123 * t246 - t250 * t162;
t51 = t124 * t250 - t126 * t246;
t122 = t163 * t247 - t251 * t209;
t276 = qJD(4) - t319 / 0.2e1;
t80 = -t197 * t334 - t198 * t335 + t204 * t251 - t247 * t205;
t166 = t211 * t246 + t242 * t347;
t261 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t69 = -pkin(4) * t316 - t80;
t17 = -pkin(4) * t308 - t27;
t260 = -t377 / 0.2e1 + t418;
t258 = t180 * mrSges(4,1) + t64 * mrSges(5,1) + Ifges(5,6) * t186 + Ifges(5,5) * t187 + t437 - (Ifges(4,4) * t248 + Ifges(4,2) * t252) * t339 / 0.2e1 - t139 * mrSges(4,3) - t65 * mrSges(5,2) + (-t282 / 0.2e1 + t382) * Ifges(5,3);
t230 = Ifges(4,5) * t252 * t315;
t229 = Ifges(5,3) * t308;
t225 = -pkin(4) - t284;
t202 = t302 * t339;
t201 = -mrSges(4,2) * t235 + mrSges(4,3) * t319;
t199 = t275 * t247;
t193 = -pkin(10) * t351 + t357;
t192 = (mrSges(4,1) * t248 + mrSges(4,2) * t252) * t315;
t172 = -t357 + (pkin(10) * t246 + pkin(5)) * t251;
t171 = -qJ(6) * t251 + t194;
t167 = t211 * t250 - t324;
t158 = mrSges(5,2) * t282 + t186 * mrSges(5,3);
t144 = t174 * t247 - t251 * t313;
t132 = -mrSges(5,2) * t308 - mrSges(5,3) * t156;
t116 = t310 + (qJD(3) * t272 + t263) * t243;
t115 = t309 + (qJD(3) * t415 + t262) * t243;
t97 = -qJD(5) * t324 + t164 * t246 + t211 * t332 - t250 * t316;
t96 = -qJD(5) * t166 + t164 * t250 + t246 * t316;
t95 = mrSges(5,1) * t156 + mrSges(5,2) * t155;
t85 = pkin(5) * t150 + qJ(6) * t149;
t82 = Ifges(5,5) * t308 - t364 + t366;
t81 = -t156 * Ifges(5,2) + Ifges(5,6) * t308 + t365;
t55 = pkin(5) * t166 - qJ(6) * t167 + t125;
t45 = qJD(4) * t123 + t115 * t247 - t251 * t312;
t44 = -qJD(4) * t122 + t115 * t251 + t247 * t312;
t43 = -pkin(5) * t210 - t51;
t42 = qJ(6) * t210 + t421;
t32 = -pkin(5) * t187 - t35;
t31 = qJ(6) * t187 + t36;
t29 = mrSges(7,1) * t78 - mrSges(7,3) * t77;
t21 = Ifges(6,4) * t77 - Ifges(6,2) * t78 + Ifges(6,6) * t156;
t18 = Ifges(7,5) * t77 + Ifges(7,6) * t156 + Ifges(7,3) * t78;
t12 = pkin(5) * t97 - qJ(6) * t96 - qJD(6) * t167 + t69;
t11 = qJD(5) * t76 - t250 * t116 + t246 * t44;
t10 = -qJD(5) * t75 + t116 * t246 + t250 * t44;
t5 = pkin(5) * t78 - qJ(6) * t77 - qJD(6) * t150 + t17;
t6 = [t115 * t201 + t123 * t132 + t44 * t158 + t162 * t95 + t209 * t192 + t380 * t76 + t381 * t75 + (-mrSges(3,1) * t249 - mrSges(3,2) * t253) * qJD(2) ^ 2 * t243 - t343 * t116 + t344 * t11 + t345 * t10 - t328 * t45 + (t29 + t360) * t122 + (t202 * t354 + (t162 * t252 - t163 * t248) * qJD(3) * mrSges(4,3)) * t339 + m(4) * (t115 * t139 - t116 * t138 + t363 + t163 * t93 + (qJD(1) * t209 + t180) * t312) + m(5) * (t116 * t120 - t122 * t27 + t123 * t26 + t44 * t65 - t45 * t64 + t363) + m(7) * (t1 * t76 + t10 * t14 + t11 * t13 + t122 * t5 + t2 * t75 + t28 * t45) + m(6) * (t10 * t25 - t11 * t24 + t122 * t17 + t3 * t76 - t4 * t75 + t45 * t53); t433 * t106 + t434 * t105 + t435 * t107 + (t1 * t42 + t2 * t43 + t5 * t55 + (t12 - t144) * t28 + t436 * t14 + t435 * t13) * m(7) + t436 * t104 + (t94 * mrSges(4,3) * t248 - pkin(2) * t192 + (-qJD(2) * t289 / 0.2e1 + t93 * mrSges(4,3) - t229 / 0.2e1 + t412) * t252 + (-m(4) * t180 - t202 + (-m(4) * pkin(2) + t302) * t339) * t321 + (((t338 / 0.2e1 + t391) * Ifges(4,5) + (-qJD(2) * t213 - t138) * mrSges(4,3) + t314) * t252 + (t437 + (-Ifges(4,6) * t244 - t214 * mrSges(4,3) + Ifges(5,5) * t392 + Ifges(5,6) * t394 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(5,3)) * t355) * qJD(2) + t258) * t248) * qJD(3) + (-0.3e1 / 0.2e1 * t248 ^ 2 + 0.3e1 / 0.2e1 * t252 ^ 2) * Ifges(4,4) * t315) * t242 + t328 * t144 + (t165 * t429 + t431 * t96) * t395 + (t167 * t431 + t210 * t429) * t401 + (t167 * t432 + t210 * t431) * t409 + (t165 * t431 + t432 * t96) * t402 + (-t164 * t64 - t165 * t65 - t210 * t26 - t211 * t27) * mrSges(5,3) + (t125 * t17 + t3 * t421 + t4 * t51 + (-t144 + t69) * t53 + t434 * t25 + t433 * t24) * m(6) + t421 * t48 + (t230 / 0.2e1 + t307) * t244 + t424 * t158 + (t140 * t27 + t141 * t26 + t196 * t94 + t424 * t65 + (t144 + t80) * t64 + t419 * t120) * m(5) + t425 * t96 / 0.2e1 + t426 * t167 / 0.2e1 + t427 * t210 / 0.2e1 + (t60 + t59) * t165 / 0.2e1 + (t1 * t210 + t14 * t165 - t167 * t5 - t28 * t96) * mrSges(7,3) + (Ifges(7,5) * t96 + Ifges(7,6) * t165) * t404 + (Ifges(7,5) * t167 + Ifges(7,6) * t210) * t407 + (-t138 * t419 + t139 * t420 - t213 * t94 + t214 * t93) * m(4) + t420 * t201 - t343 * t419 - t156 * (Ifges(5,4) * t211 - Ifges(5,2) * t210) / 0.2e1 + t186 * (Ifges(5,4) * t164 - Ifges(5,2) * t165) / 0.2e1 + t81 * t394 + t82 * t392 + t289 * t382 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) + t5 * mrSges(7,1) + t17 * mrSges(6,1) + t18 / 0.2e1 - t21 / 0.2e1 + Ifges(7,3) * t407 - Ifges(6,2) * t408 + t430 * t409 - t447 * t401) * t166 + (-Ifges(6,2) * t405 + Ifges(7,3) * t404 - t395 * t447 + t430 * t402 + t442) * t97 + (Ifges(6,4) * t96 + Ifges(6,6) * t165) * t405 + (Ifges(6,4) * t167 + Ifges(6,6) * t210) * t408 + t165 * t441 + t2 * (-mrSges(7,1) * t210 + mrSges(7,2) * t167) + t4 * (mrSges(6,1) * t210 - mrSges(6,3) * t167) + t94 * (mrSges(5,1) * t210 + mrSges(5,2) * t211) + t196 * t95 + t24 * (mrSges(6,1) * t165 - mrSges(6,3) * t96) + t120 * (mrSges(5,1) * t165 + mrSges(5,2) * t164) + t13 * (-mrSges(7,1) * t165 + mrSges(7,2) * t96) + t164 * t113 / 0.2e1 + t80 * t159 + t140 * t131 + t141 * t132 + t125 * t30 + t155 * (Ifges(5,1) * t211 - Ifges(5,4) * t210) / 0.2e1 + t187 * (Ifges(5,1) * t164 - Ifges(5,4) * t165) / 0.2e1 + (-t165 * t25 + t167 * t17 - t210 * t3 + t53 * t96) * mrSges(6,2) + t43 * t47 + t42 * t49 + t51 * t46 + t55 * t29 + t12 * t86 + t69 * t87; (-m(5) * t94 - t95) * pkin(3) + t444 * t104 + t445 * t86 + (t1 * t171 + t13 * t446 + t14 * t444 + t172 * t2 + t199 * t5 + t28 * t445) * m(7) + t446 * t107 + t343 * t139 + m(6) * (t135 * t25 + t136 * t24 + t193 * t4 + t194 * t3) + (t26 * mrSges(5,3) - t20 / 0.2e1 + t81 / 0.2e1 - t94 * mrSges(5,1) + t365 / 0.2e1 - t19 / 0.2e1 - t325 * t78 + t327 * t77 + (m(5) * t26 + t132) * pkin(10) + (-Ifges(5,2) / 0.2e1 - t326) * t156 + t261) * t251 + (-Ifges(6,2) * t404 + Ifges(7,3) * t405 - t396 * t447 + t403 * t430 - t442) * t178 + (-mrSges(6,2) * t53 - t13 * mrSges(7,2) + t24 * mrSges(6,3) + mrSges(7,3) * t28 + Ifges(6,4) * t404 + Ifges(7,5) * t405 + t396 * t431 + t403 * t432 - t425 / 0.2e1) * t179 - m(6) * (t24 * t37 + t25 * t38 + t53 * t90) - m(5) * (t101 * t64 + t102 * t65 + t120 * t139) + t451 * t106 + t452 * t105 + ((-Ifges(5,6) * qJD(4) / 0.2e1 + (-m(5) * t65 - t158) * pkin(10) - t443) * t247 + (t255 + Ifges(5,5) * t382 + (-m(5) * t64 + m(6) * t53 + t359) * pkin(10) - t260) * t251) * qJD(4) + ((qJD(3) * (Ifges(5,5) * t247 + Ifges(5,6) * t251) / 0.2e1 + Ifges(4,4) * t438 + (-qJD(3) + t391) * Ifges(4,6) - t258) * t248 + (t138 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t320 + (-Ifges(5,5) * t276 + t260) * t251 + (t276 * Ifges(5,6) + t443) * t247 - t314 - t411) * t252) * t339 + t199 * t29 - t138 * t201 + t193 * t46 + t194 * t48 + t171 * t49 + t172 * t47 - t102 * t158 - t101 * t159 + (t366 / 0.2e1 - t364 / 0.2e1 + t94 * mrSges(5,2) + t82 / 0.2e1 + t17 * t300 + t5 * t298 - t27 * mrSges(5,3) + t18 * t389 + t21 * t390 + t286 * t407 + t293 * t408 + (-t246 * t3 - t250 * t4) * mrSges(6,3) + (-t1 * t246 + t2 * t250) * mrSges(7,2) + (-m(5) * t27 + m(6) * t17 + t360) * pkin(10) + (t285 * t405 + t292 * t404 + t28 * t299 + t53 * t301 + t61 * t388 + (t24 * t246 - t25 * t250) * mrSges(6,3) + (-t13 * t246 - t14 * t250) * mrSges(7,2) + t416 * t403 + t417 * t396 + t425 * t390) * qJD(5) + t439 * t409 + t440 * t401 + (qJD(5) * t58 + t426) * t387) * t247 + t307 + t230 - t90 * t87; -t359 * t65 + t255 * qJD(5) + (-t13 * t32 - t14 * t31 + t225 * t5 + t28 * t423) * m(7) + t423 * t86 + t426 * t389 - t5 * t299 - t17 * t301 + t413 * mrSges(7,2) + t414 * mrSges(6,3) + ((-m(6) * t279 + m(7) * t281 - t246 * t345 + t250 * t344) * qJD(5) + m(6) * t414 + m(7) * t413 + t246 * t381 + t250 * t380) * pkin(11) + t416 * t409 + t417 * t401 + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t187 + t270 / 0.2e1 - t255 + t418) * t186 + t21 * t387 + t18 * t388 + (-t410 - t269 / 0.2e1) * t187 + t225 * t29 - t64 * t158 + t285 * t407 + t292 * t408 - t412 + t229 + (-pkin(4) * t17 - t24 * t35 - t25 * t36 - t53 * t65) * m(6) - pkin(4) * t30 - t31 * t104 - t36 * t105 - t35 * t106 - t32 * t107; -t261 + (t13 * t149 + t14 * t150) * mrSges(7,2) + (-t344 + t378) * t25 + (-t345 - t379) * t24 - t53 * (mrSges(6,1) * t150 - mrSges(6,2) * t149) + t61 * t402 + (Ifges(7,3) * t150 - t371) * t405 - t28 * (mrSges(7,1) * t150 + mrSges(7,3) * t149) - pkin(5) * t47 + qJ(6) * t49 - t85 * t86 + qJD(6) * t104 + (-t149 * t431 - t150 * t447) * t396 + (-pkin(5) * t2 + qJ(6) * t1 - t13 * t25 + t14 * t422 - t28 * t85) * m(7) + (-Ifges(6,2) * t150 - t148 + t425) * t404 + (-t149 * t432 + t147 - t375 + t58) * t403 + t427; -t182 * t104 + t150 * t86 + 0.2e1 * (t2 / 0.2e1 + t14 * t396 + t28 * t402) * m(7) + t47;];
tauc  = t6(:);

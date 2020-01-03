% Calculate vector of inverse dynamics joint torques for
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR16_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:44
% EndTime: 2019-12-31 20:45:21
% DurationCPUTime: 19.74s
% Computational Cost: add. (6852->723), mult. (16691->975), div. (0->0), fcn. (12246->10), ass. (0->335)
t232 = cos(qJ(4));
t391 = pkin(2) + pkin(8);
t327 = qJD(4) * t391;
t229 = sin(qJ(2));
t225 = sin(pkin(5));
t335 = qJD(1) * t225;
t310 = t229 * t335;
t201 = pkin(2) * t310;
t233 = cos(qJ(2));
t260 = pkin(8) * t229 - qJ(3) * t233;
t125 = t260 * t335 + t201;
t309 = t233 * t335;
t204 = pkin(7) * t309;
t226 = cos(pkin(5));
t334 = qJD(1) * t226;
t320 = pkin(1) * t334;
t161 = t229 * t320 + t204;
t127 = pkin(3) * t309 + t161;
t228 = sin(qJ(4));
t64 = t232 * t125 + t228 * t127;
t448 = -t232 * t327 - t64;
t227 = sin(qJ(5));
t231 = cos(qJ(5));
t271 = t227 * mrSges(6,1) + t231 * mrSges(6,2);
t392 = -m(6) - m(5);
t445 = -mrSges(3,1) + mrSges(4,2);
t418 = pkin(8) * t392 - mrSges(5,3) + t445;
t447 = t271 - t418;
t206 = t233 * t320;
t446 = qJD(3) - t206;
t278 = pkin(4) * t232 + pkin(9) * t228;
t390 = pkin(3) + pkin(7);
t444 = qJD(4) * t278 - (-t278 - t390) * t310 + t446;
t443 = pkin(9) * t309 - t448;
t211 = qJD(2) + t334;
t140 = -t211 * t228 - t232 * t309;
t137 = qJD(5) - t140;
t187 = qJD(4) + t310;
t245 = -t211 * t232 + t228 * t309;
t87 = t187 * t231 + t227 * t245;
t88 = t187 * t227 - t231 * t245;
t35 = t88 * Ifges(6,5) + t87 * Ifges(6,6) + t137 * Ifges(6,3);
t85 = -t211 * t391 + t310 * t390 + t446;
t297 = -qJ(3) * t229 - pkin(1);
t98 = (-t233 * t391 + t297) * t335;
t42 = t228 * t85 + t232 * t98;
t39 = pkin(9) * t187 + t42;
t193 = t211 * qJ(3);
t93 = t193 + t127;
t43 = -pkin(4) * t140 + pkin(9) * t245 + t93;
t15 = t227 * t43 + t231 * t39;
t433 = t15 * mrSges(6,2);
t14 = -t227 * t39 + t231 * t43;
t434 = t14 * mrSges(6,1);
t368 = Ifges(5,4) * t245;
t439 = t187 * Ifges(5,6);
t441 = t140 * Ifges(5,2);
t67 = -t368 + t439 + t441;
t442 = -t433 + t434 - t42 * mrSges(5,3) + t35 / 0.2e1 - t67 / 0.2e1;
t323 = m(4) - t392;
t135 = Ifges(5,4) * t140;
t440 = t187 * Ifges(5,5);
t349 = t225 * t233;
t272 = -mrSges(6,1) * t231 + mrSges(6,2) * t227;
t248 = m(6) * pkin(4) - t272;
t274 = mrSges(5,1) * t228 + mrSges(5,2) * t232;
t438 = -t248 * t228 - t274;
t351 = t225 * t229;
t308 = qJD(2) * t351;
t164 = qJD(1) * t308 - qJDD(1) * t349;
t437 = pkin(2) * t323 + t447;
t322 = qJDD(1) * t226;
t209 = qJDD(2) + t322;
t381 = pkin(1) * t226;
t319 = qJD(2) * t381;
t288 = qJD(1) * t319;
t315 = pkin(1) * t322;
t89 = -pkin(7) * t164 + t229 * t315 + t233 * t288;
t72 = -t209 * qJ(3) - t211 * qJD(3) - t89;
t50 = -pkin(3) * t164 - t72;
t73 = qJD(4) * t140 + t164 * t228 + t209 * t232;
t74 = qJD(4) * t245 + t164 * t232 - t209 * t228;
t13 = -pkin(4) * t74 - pkin(9) * t73 + t50;
t328 = qJD(4) * t232;
t329 = qJD(4) * t228;
t331 = qJD(2) * t233;
t165 = (qJD(1) * t331 + qJDD(1) * t229) * t225;
t212 = pkin(7) * t351;
t90 = -qJD(2) * t204 - qJDD(1) * t212 - t229 * t288 + t233 * t315;
t241 = qJDD(3) - t90;
t49 = pkin(3) * t165 - t209 * t391 + t241;
t330 = qJD(3) * t229;
t239 = -qJ(3) * t165 + (-pkin(1) * qJDD(1) - qJD(1) * t330) * t225;
t53 = t164 * t391 + t239;
t11 = t228 * t49 + t232 * t53 + t85 * t328 - t329 * t98;
t150 = qJDD(4) + t165;
t8 = pkin(9) * t150 + t11;
t1 = qJD(5) * t14 + t13 * t227 + t231 * t8;
t2 = -qJD(5) * t15 + t13 * t231 - t227 * t8;
t436 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t384 = t150 / 0.2e1;
t397 = t74 / 0.2e1;
t398 = t73 / 0.2e1;
t408 = Ifges(5,1) * t398 + Ifges(5,4) * t397 + Ifges(5,5) * t384;
t32 = qJD(5) * t87 + t150 * t227 + t231 * t73;
t407 = t32 / 0.2e1;
t33 = -qJD(5) * t88 + t150 * t231 - t227 * t73;
t406 = t33 / 0.2e1;
t71 = qJDD(5) - t74;
t399 = t71 / 0.2e1;
t374 = mrSges(4,3) - mrSges(3,2);
t432 = -Ifges(4,4) + Ifges(3,5);
t431 = Ifges(4,5) - Ifges(3,6);
t10 = -mrSges(6,1) * t33 + mrSges(6,2) * t32;
t46 = mrSges(5,1) * t150 - mrSges(5,3) * t73;
t430 = -t10 + t46;
t295 = -pkin(9) * t232 + qJ(3);
t180 = pkin(4) * t228 + t295;
t345 = t228 * t391;
t139 = t227 * t180 - t231 * t345;
t429 = -qJD(5) * t139 + t227 * t443 + t231 * t444;
t138 = t231 * t180 + t227 * t345;
t428 = qJD(5) * t138 + t227 * t444 - t231 * t443;
t44 = -mrSges(6,1) * t87 + mrSges(6,2) * t88;
t370 = mrSges(5,3) * t245;
t95 = mrSges(5,1) * t187 + t370;
t427 = t44 - t95;
t426 = mrSges(5,1) + t248;
t316 = m(6) * pkin(9) + mrSges(6,3);
t294 = mrSges(5,2) - t316;
t291 = mrSges(4,1) * t309;
t156 = -mrSges(4,3) * t211 - t291;
t83 = -mrSges(5,1) * t140 - mrSges(5,2) * t245;
t425 = -t83 + t156;
t346 = t228 * t229;
t133 = (-t227 * t346 + t231 * t233) * t335;
t306 = t227 * t329;
t324 = qJD(5) * t232;
t424 = t231 * t324 + t133 - t306;
t343 = t229 * t231;
t134 = (t227 * t233 + t228 * t343) * t335;
t423 = t227 * t324 + t231 * t329 + t134;
t290 = mrSges(3,3) * t310;
t292 = mrSges(4,1) * t310;
t422 = t211 * t445 + t290 + t292;
t16 = mrSges(6,1) * t71 - mrSges(6,3) * t32;
t17 = -mrSges(6,2) * t71 + mrSges(6,3) * t33;
t421 = -t227 * t16 + t231 * t17;
t420 = t1 * t231 - t2 * t227;
t369 = Ifges(3,4) * t229;
t419 = -t229 * (Ifges(3,1) * t233 - t369) / 0.2e1 + pkin(1) * (mrSges(3,1) * t229 + mrSges(3,2) * t233);
t416 = qJ(3) * t323 + t374;
t12 = -qJD(4) * t42 - t228 * t53 + t232 * t49;
t415 = t12 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,5) * t73 + Ifges(5,6) * t74 + Ifges(5,3) * t150;
t221 = t229 * t381;
t128 = (t349 * t390 + t221) * qJD(2);
t380 = pkin(1) * t233;
t311 = -pkin(2) - t380;
t100 = pkin(3) * t351 + t212 + (-pkin(8) + t311) * t226;
t337 = pkin(2) * t349 + qJ(3) * t351;
t382 = pkin(1) * t225;
t144 = -t337 - t382;
t215 = pkin(8) * t349;
t119 = t144 - t215;
t354 = t228 * t100 + t232 * t119;
t203 = pkin(2) * t308;
t97 = t203 + (qJD(2) * t260 - t330) * t225;
t27 = -qJD(4) * t354 + t128 * t232 - t228 * t97;
t413 = -m(6) * t295 + t232 * mrSges(6,3) - t374 + t438 + (-m(4) - m(5)) * qJ(3);
t5 = Ifges(6,5) * t32 + Ifges(6,6) * t33 + Ifges(6,3) * t71;
t412 = t5 / 0.2e1;
t6 = t32 * Ifges(6,4) + t33 * Ifges(6,2) + t71 * Ifges(6,6);
t411 = t6 / 0.2e1;
t410 = Ifges(6,1) * t407 + Ifges(6,4) * t406 + Ifges(6,5) * t399;
t409 = -t73 * Ifges(5,4) / 0.2e1 - t74 * Ifges(5,2) / 0.2e1 - t150 * Ifges(5,6) / 0.2e1;
t383 = Ifges(6,4) * t88;
t36 = Ifges(6,2) * t87 + Ifges(6,6) * t137 + t383;
t404 = -t36 / 0.2e1;
t403 = t36 / 0.2e1;
t86 = Ifges(6,4) * t87;
t37 = Ifges(6,1) * t88 + Ifges(6,5) * t137 + t86;
t402 = -t37 / 0.2e1;
t401 = t37 / 0.2e1;
t396 = -t87 / 0.2e1;
t395 = t87 / 0.2e1;
t394 = -t88 / 0.2e1;
t393 = t88 / 0.2e1;
t389 = -t137 / 0.2e1;
t388 = t137 / 0.2e1;
t385 = -t245 / 0.2e1;
t9 = -pkin(4) * t150 - t12;
t377 = t232 * t9;
t373 = Ifges(3,4) + Ifges(4,6);
t371 = mrSges(5,3) * t140;
t367 = Ifges(5,4) * t228;
t366 = Ifges(5,4) * t232;
t365 = Ifges(6,4) * t227;
t364 = Ifges(6,4) * t231;
t363 = Ifges(4,6) * t229;
t362 = Ifges(4,6) * t233;
t361 = t11 * t228;
t360 = t140 * Ifges(5,6);
t359 = t245 * Ifges(5,5);
t358 = t187 * Ifges(5,3);
t355 = t232 * mrSges(5,3);
t353 = t140 * t227;
t352 = t140 * t231;
t230 = sin(qJ(1));
t350 = t225 * t230;
t234 = cos(qJ(1));
t348 = t225 * t234;
t347 = t227 * t232;
t344 = t229 * t230;
t342 = t229 * t234;
t341 = t230 * t233;
t340 = t231 * t232;
t338 = t233 * t234;
t178 = pkin(7) * t349 + t221;
t336 = t234 * pkin(1) + pkin(7) * t350;
t332 = qJD(1) ^ 2 * t225 ^ 2;
t326 = qJD(5) * t227;
t325 = qJD(5) * t231;
t318 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t317 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t314 = t228 * t349;
t175 = -t226 * t344 + t338;
t313 = t175 * pkin(2) + t336;
t143 = -t226 * qJ(3) - t178;
t307 = t225 * t331;
t299 = -t329 / 0.2e1;
t298 = -t324 / 0.2e1;
t296 = -pkin(1) * t230 + pkin(7) * t348;
t106 = t165 * mrSges(4,1) + t209 * mrSges(4,2);
t293 = t390 * t351;
t289 = mrSges(3,3) * t309;
t287 = pkin(3) * t350 + t313;
t118 = pkin(3) * t349 - t143;
t286 = t232 * t310;
t173 = t226 * t342 + t341;
t280 = -t173 * pkin(2) + t296;
t170 = t226 * t228 + t232 * t349;
t171 = t226 * t232 - t314;
t276 = mrSges(5,1) * t170 + mrSges(5,2) * t171;
t109 = -t171 * t227 + t225 * t343;
t110 = t171 * t231 + t227 * t351;
t273 = mrSges(6,1) * t109 - mrSges(6,2) * t110;
t270 = mrSges(4,2) * t233 - mrSges(4,3) * t229;
t269 = Ifges(5,1) * t228 + t366;
t268 = Ifges(6,1) * t231 - t365;
t267 = Ifges(6,1) * t227 + t364;
t266 = Ifges(5,2) * t232 + t367;
t265 = -Ifges(6,2) * t227 + t364;
t264 = Ifges(6,2) * t231 + t365;
t263 = Ifges(5,5) * t228 + Ifges(5,6) * t232;
t262 = Ifges(6,5) * t231 - Ifges(6,6) * t227;
t261 = Ifges(6,5) * t227 + Ifges(6,6) * t231;
t55 = pkin(9) * t351 + t354;
t62 = pkin(4) * t170 - pkin(9) * t171 + t118;
t21 = t227 * t62 + t231 * t55;
t20 = -t227 * t55 + t231 * t62;
t41 = -t228 * t98 + t232 * t85;
t259 = t228 * t41 - t232 * t42;
t58 = t100 * t232 - t119 * t228;
t63 = -t125 * t228 + t127 * t232;
t160 = pkin(7) * t310 - t206;
t207 = t233 * t319;
t162 = -pkin(7) * t308 + t207;
t172 = -t226 * t338 + t344;
t115 = -t172 * t228 + t232 * t348;
t113 = t172 * t232 + t228 * t348;
t250 = t229 * (-Ifges(4,2) * t233 + t363);
t249 = t233 * (Ifges(4,3) * t229 - t362);
t26 = t100 * t328 - t119 * t329 + t228 * t128 + t232 * t97;
t222 = t226 * qJD(3);
t99 = -qJD(2) * t293 + t207 + t222;
t81 = -pkin(2) * t209 + t241;
t240 = t90 * mrSges(3,1) - t89 * mrSges(3,2) + t81 * mrSges(4,2) - t72 * mrSges(4,3);
t238 = (-t14 * t231 - t15 * t227) * qJD(5) + t420;
t237 = -qJD(4) * t259 + t12 * t232 + t361;
t200 = Ifges(3,4) * t309;
t192 = Ifges(4,1) * t209;
t191 = Ifges(3,3) * t209;
t177 = t226 * t380 - t212;
t176 = (-mrSges(3,1) * t233 + mrSges(3,2) * t229) * t225;
t174 = t226 * t341 + t342;
t163 = t178 * qJD(2);
t159 = -qJ(3) * t309 + t201;
t158 = t270 * t335;
t155 = -mrSges(3,2) * t211 + t289;
t149 = Ifges(4,4) * t165;
t148 = Ifges(3,5) * t165;
t147 = Ifges(4,5) * t164;
t146 = Ifges(3,6) * t164;
t145 = t226 * t311 + t212;
t136 = -t162 - t222;
t132 = t211 * t227 - t231 * t286;
t131 = t211 * t231 + t227 * t286;
t130 = (-pkin(2) * t233 + t297) * t335;
t129 = t203 + (-qJ(3) * t331 - t330) * t225;
t126 = -qJD(1) * t293 + t206;
t124 = -t193 - t161;
t123 = t211 * Ifges(4,4) + (-Ifges(4,2) * t229 - t362) * t335;
t122 = t211 * Ifges(4,5) + (-t233 * Ifges(4,3) - t363) * t335;
t121 = Ifges(3,1) * t310 + t211 * Ifges(3,5) + t200;
t120 = t211 * Ifges(3,6) + (t233 * Ifges(3,2) + t369) * t335;
t117 = -pkin(2) * t211 + qJD(3) + t160;
t112 = t174 * t228 + t232 * t350;
t111 = -t174 * t232 + t228 * t350;
t108 = -qJD(4) * t314 + t226 * t328 - t232 * t308;
t107 = -qJD(4) * t170 + t228 * t308;
t105 = mrSges(4,1) * t164 - mrSges(4,3) * t209;
t94 = -mrSges(5,2) * t187 + t371;
t84 = -pkin(4) * t245 - pkin(9) * t140;
t79 = t112 * t231 + t175 * t227;
t78 = -t112 * t227 + t175 * t231;
t75 = pkin(2) * t164 + t239;
t68 = -Ifges(5,1) * t245 + t135 + t440;
t66 = t358 - t359 + t360;
t61 = mrSges(6,1) * t137 - mrSges(6,3) * t88;
t60 = -mrSges(6,2) * t137 + mrSges(6,3) * t87;
t56 = -pkin(4) * t309 - t63;
t54 = -pkin(4) * t351 - t58;
t52 = qJD(5) * t109 + t107 * t231 + t227 * t307;
t51 = -qJD(5) * t110 - t107 * t227 + t231 * t307;
t47 = -mrSges(5,2) * t150 + mrSges(5,3) * t74;
t40 = pkin(4) * t108 - pkin(9) * t107 + t99;
t38 = -pkin(4) * t187 - t41;
t34 = -mrSges(5,1) * t74 + mrSges(5,2) * t73;
t23 = t227 * t84 + t231 * t41;
t22 = -t227 * t41 + t231 * t84;
t19 = -pkin(4) * t307 - t27;
t18 = pkin(9) * t307 + t26;
t4 = -qJD(5) * t21 - t18 * t227 + t231 * t40;
t3 = qJD(5) * t20 + t18 * t231 + t227 * t40;
t7 = [(t1 * t109 - t110 * t2 - t14 * t52 + t15 * t51) * mrSges(6,3) + (t135 / 0.2e1 - t41 * mrSges(5,3) + t440 / 0.2e1 + Ifges(5,1) * t385 + t68 / 0.2e1 + t93 * mrSges(5,2)) * t107 + (Ifges(6,1) * t52 + Ifges(6,4) * t51) * t393 + (Ifges(6,1) * t110 + Ifges(6,4) * t109) * t407 + (Ifges(6,5) * t52 + Ifges(6,6) * t51) * t388 + (Ifges(6,5) * t110 + Ifges(6,6) * t109) * t399 + m(6) * (t1 * t21 + t14 * t4 + t15 * t3 + t19 * t38 + t2 * t20 + t54 * t9) + m(4) * (t117 * t163 + t124 * t136 + t129 * t130 + t143 * t72 + t144 * t75 + t145 * t81) + (Ifges(6,4) * t52 + Ifges(6,2) * t51) * t395 + (Ifges(6,4) * t110 + Ifges(6,2) * t109) * t406 + (-m(4) * t280 - m(3) * t296 + mrSges(2,1) * t230 + mrSges(2,2) * t234 + t416 * t172 + t294 * t113 - t426 * t115 + t447 * t173 - t392 * (-pkin(3) * t348 - t280)) * g(1) + (-m(4) * t313 - m(3) * t336 - mrSges(2,1) * t234 + mrSges(2,2) * t230 - m(6) * (pkin(4) * t112 + t287) - t79 * mrSges(6,1) - t78 * mrSges(6,2) - m(5) * t287 - t112 * mrSges(5,1) - t416 * t174 + t294 * t111 + t418 * t175) * g(2) + m(3) * (t160 * t163 + t161 * t162 + t177 * t90 + t178 * t89) + (-t12 * mrSges(5,3) + 0.2e1 * t408) * t171 + (-mrSges(5,3) * t11 - Ifges(5,4) * t398 + Ifges(6,5) * t407 - Ifges(5,2) * t397 - Ifges(5,6) * t384 + Ifges(6,6) * t406 + Ifges(6,3) * t399 + t409 + t412 + t436) * t170 + t50 * t276 + (t148 / 0.2e1 - t146 / 0.2e1 + t191 / 0.2e1 + t192 / 0.2e1 - t149 / 0.2e1 + t147 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t209 + t318 * t165 + t317 * t164 + t240) * t226 + t422 * t163 + (-t441 / 0.2e1 + Ifges(6,5) * t393 - t439 / 0.2e1 - Ifges(5,4) * t385 + Ifges(6,3) * t388 + Ifges(6,6) * t395 + t93 * mrSges(5,1) + t442) * t108 + m(5) * (t11 * t354 + t118 * t50 + t12 * t58 + t26 * t42 + t27 * t41 + t93 * t99) + t354 * t47 + ((-mrSges(3,1) * t164 - mrSges(3,2) * t165 + (m(3) * t382 - t176) * qJDD(1)) * pkin(1) + (-t72 * mrSges(4,1) + t75 * mrSges(4,2) + t89 * mrSges(3,3) - t431 * t209 + t373 * t165 + (-Ifges(3,2) - Ifges(4,3)) * t164) * t233 + (t81 * mrSges(4,1) - t90 * mrSges(3,3) - t75 * mrSges(4,3) + t432 * t209 + (Ifges(4,2) + Ifges(3,1)) * t165 - t373 * t164 + t415) * t229 + ((-t120 / 0.2e1 + t122 / 0.2e1 - t130 * mrSges(4,2) + t124 * mrSges(4,1) - t161 * mrSges(3,3) + t317 * t211) * t229 + (t121 / 0.2e1 - t123 / 0.2e1 + t66 / 0.2e1 - t130 * mrSges(4,3) + t160 * mrSges(3,3) + t117 * mrSges(4,1) + t41 * mrSges(5,1) - t42 * mrSges(5,2) + t360 / 0.2e1 - t359 / 0.2e1 + t358 / 0.2e1 + t318 * t211) * t233 + (t233 * (Ifges(3,4) * t233 - Ifges(3,2) * t229) / 0.2e1 - t249 / 0.2e1 - t250 / 0.2e1 - t419) * t335) * qJD(2) + (mrSges(3,3) + mrSges(4,1)) * (-g(1) * t234 - g(2) * t230)) * t225 - t9 * t273 + Ifges(2,3) * qJDD(1) + t177 * (mrSges(3,1) * t209 - mrSges(3,3) * t165) + t178 * (-mrSges(3,2) * t209 - mrSges(3,3) * t164) + t144 * (-mrSges(4,2) * t164 - mrSges(4,3) * t165) + t129 * t158 + t162 * t155 + t136 * t156 + t52 * t401 + t51 * t403 + t110 * t410 + t109 * t411 + t20 * t16 + t21 * t17 + t19 * t44 + t38 * (-mrSges(6,1) * t51 + mrSges(6,2) * t52) + t54 * t10 + t58 * t46 + t3 * t60 + t4 * t61 + t26 * t94 + t27 * t95 + t99 * t83 + t118 * t34 + t143 * t105 + t145 * t106; -((t229 * t431 + t233 * t432) * t211 + (-Ifges(3,2) * t310 + t121 + t200 + t66) * t233 + (t228 * t68 + t232 * t67 + t122) * t229 + t187 * (Ifges(5,3) * t233 + t229 * t263) - t245 * (Ifges(5,5) * t233 + t229 * t269) + t140 * (Ifges(5,6) * t233 + t229 * t266)) * t335 / 0.2e1 - (t140 * t266 + t187 * t263 - t245 * t269) * qJD(4) / 0.2e1 + (-t38 * t56 + t1 * t139 + t2 * t138 - (t329 * t38 - t377) * t391 + t428 * t15 + t429 * t14) * m(6) + (t50 * qJ(3) - t391 * t237 - t41 * t63 - t42 * t64 + (qJD(3) - t126) * t93) * m(5) + (t233 * t123 + (t232 * t35 + t120) * t229) * t335 / 0.2e1 + t187 * t93 * (mrSges(5,1) * t232 - mrSges(5,2) * t228) - t427 * t228 * t327 - t430 * t232 * t391 - t6 * t347 / 0.2e1 + (-pkin(2) * t81 - qJ(3) * t72 - qJD(3) * t124 - t130 * t159) * m(4) + t148 - t149 - t146 + t147 + t240 - t124 * t292 + (-t267 * t324 + (Ifges(6,5) * t232 - t228 * t268) * qJD(4)) * t393 + (Ifges(6,1) * t134 + Ifges(6,4) * t133 - Ifges(6,5) * t286) * t394 + (-t264 * t324 + (Ifges(6,6) * t232 - t228 * t265) * qJD(4)) * t395 + t227 * t37 * t298 - t117 * t291 + (t34 - t105) * qJ(3) + t50 * t274 + t419 * t332 + (t298 * t36 + t299 * t37) * t231 + (t174 * t437 + t175 * t413) * g(1) + (t172 * t437 + t173 * t413) * g(2) + (-m(4) * t337 + t176 + t392 * (t215 + t337) + (t270 + (-t271 - mrSges(5,3)) * t233 + (t232 * t316 + t438) * t229) * t225) * g(3) + (-m(4) * t117 + t290 - t422) * t161 + (-t1 * t228 - t15 * t286 - t38 * t423) * mrSges(6,2) + (t14 * t286 + t2 * t228 + t38 * t424) * mrSges(6,1) + (-t1 * t347 + t14 * t423 - t15 * t424 - t2 * t340) * mrSges(6,3) + t448 * t94 + t442 * t328 + (Ifges(5,5) * t232 - Ifges(5,6) * t228) * t384 + (-t261 * t324 + (Ifges(6,3) * t232 - t228 * t262) * qJD(4)) * t388 + (Ifges(6,5) * t134 + Ifges(6,6) * t133 - Ifges(6,3) * t286) * t389 + t271 * t377 - t47 * t345 + t428 * t60 + t429 * t61 + (t250 + t249) * t332 / 0.2e1 - t12 * t355 + t191 + t192 + (-m(4) * t124 + t155 - t156 - t289) * t160 - t159 * t158 + (Ifges(6,4) * t134 + Ifges(6,2) * t133 - Ifges(6,6) * t286) * t396 + (-Ifges(5,2) * t228 + t366) * t397 + (Ifges(5,1) * t232 - t367) * t398 + (Ifges(6,3) * t228 + t232 * t262) * t399 + (t329 * t41 - t361) * mrSges(5,3) + t68 * t299 + t134 * t402 + t306 * t403 + t133 * t404 + (Ifges(6,6) * t228 + t232 * t265) * t406 + (Ifges(6,5) * t228 + t232 * t268) * t407 + t232 * t408 + t228 * t409 + t340 * t410 + t228 * t412 - t425 * qJD(3) + (-t41 * (mrSges(5,1) * t233 - mrSges(5,3) * t346) - t42 * (-mrSges(5,2) * t233 + t229 * t355) - t130 * (-mrSges(4,2) * t229 - mrSges(4,3) * t233)) * t335 - t56 * t44 - t63 * t95 - pkin(2) * t106 - t126 * t83 + t138 * t16 + t139 * t17; -t131 * t61 - t132 * t60 + t425 * t211 + t158 * t310 + (t94 * t310 + (-t227 * t61 + t231 * t60 + t94) * qJD(4) + t430) * t232 + (t47 + (-t227 * t60 - t231 * t61) * qJD(5) + t187 * t427 + t421) * t228 + t106 + (-t174 * g(1) - t172 * g(2) + g(3) * t349) * t323 + (-t131 * t14 - t132 * t15 + (-t9 + (-t14 * t227 + t15 * t231) * qJD(4)) * t232 + (t187 * t38 + t238) * t228) * m(6) + (-t211 * t93 - t259 * t310 + t237) * m(5) + (t124 * t211 + t130 * t310 + t81) * m(4); (-Ifges(6,3) * t245 + t140 * t262) * t389 - (Ifges(5,2) * t245 + t135 + t68) * t140 / 0.2e1 - t187 * (Ifges(5,5) * t140 + Ifges(5,6) * t245) / 0.2e1 + (-Ifges(6,6) * t245 + t140 * t265) * t396 - t93 * (-mrSges(5,1) * t245 + mrSges(5,2) * t140) + (Ifges(5,1) * t140 + t35 + t368) * t245 / 0.2e1 + t245 * t434 - t245 * t433 + (-Ifges(6,5) * t245 + t140 * t268) * t394 + t137 * t38 * t271 + t415 + (t170 * t248 - t171 * t316 + t276) * g(3) + (-pkin(4) * t9 - t14 * t22 - t15 * t23) * m(6) + (-t113 * t426 - t115 * t294) * g(2) + ((-t326 + t353) * t15 + (-t325 + t352) * t14 + t420) * mrSges(6,3) + (m(6) * t238 - t325 * t61 - t326 * t60 + t421) * pkin(9) + t67 * t385 + (-m(6) * t38 - t370 - t427) * t42 + (t137 * t262 + t265 * t87 + t268 * t88) * qJD(5) / 0.2e1 + t9 * t272 + t261 * t399 + t325 * t401 + t352 * t402 + t353 * t403 + t326 * t404 + t264 * t406 + t267 * t407 + t227 * t410 + t231 * t411 + (t371 - t94) * t41 - pkin(4) * t10 + (t111 * t426 + t112 * t294) * g(1) - t23 * t60 - t22 * t61; -t38 * (mrSges(6,1) * t88 + mrSges(6,2) * t87) + (Ifges(6,1) * t87 - t383) * t394 + t36 * t393 + (Ifges(6,5) * t87 - Ifges(6,6) * t88) * t389 - t14 * t60 + t15 * t61 - g(1) * (mrSges(6,1) * t78 - mrSges(6,2) * t79) - g(2) * ((t115 * t227 + t173 * t231) * mrSges(6,1) + (t115 * t231 - t173 * t227) * mrSges(6,2)) - g(3) * t273 + (t14 * t87 + t15 * t88) * mrSges(6,3) + t5 + (-Ifges(6,2) * t88 + t37 + t86) * t396 + t436;];
tau = t7;

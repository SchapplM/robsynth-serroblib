% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:10:02
% EndTime: 2019-03-09 08:10:12
% DurationCPUTime: 5.80s
% Computational Cost: add. (14719->426), mult. (28012->578), div. (0->0), fcn. (31759->8), ass. (0->234)
t251 = sin(pkin(9));
t254 = sin(qJ(2));
t342 = cos(pkin(9));
t377 = cos(qJ(2));
t229 = t251 * t377 + t254 * t342;
t250 = sin(pkin(10));
t252 = cos(pkin(10));
t253 = sin(qJ(6));
t376 = cos(qJ(6));
t279 = t250 * t376 + t253 * t252;
t270 = t279 * t229;
t278 = -t253 * t250 + t252 * t376;
t408 = t278 * t229;
t425 = m(7) * (t278 * t270 - t279 * t408);
t235 = (-qJ(3) - pkin(7)) * t254;
t312 = t377 * pkin(7);
t236 = qJ(3) * t377 + t312;
t198 = -t342 * t235 + t236 * t251;
t302 = -pkin(5) * t252 - pkin(4);
t103 = t229 * t302 - t198;
t149 = pkin(4) * t229 + t198;
t392 = -m(7) / 0.2e1;
t394 = -m(6) / 0.2e1;
t356 = t270 * mrSges(7,2);
t358 = t408 * mrSges(7,1);
t407 = t358 / 0.2e1 - t356 / 0.2e1;
t424 = -t103 * t392 + t149 * t394 - t407;
t225 = t251 * t254 - t342 * t377;
t245 = -pkin(2) * t377 - pkin(1);
t272 = -t229 * qJ(4) + t245;
t177 = t225 * pkin(3) + t272;
t423 = m(5) * t177 - mrSges(5,2) * t225 - mrSges(5,3) * t229;
t391 = m(7) / 0.2e1;
t421 = -mrSges(4,2) + mrSges(5,3);
t157 = t278 * t225;
t160 = t279 * t225;
t313 = (t157 * t279 - t160 * t278) * t391;
t247 = t254 * pkin(2);
t297 = qJ(4) * t225 + t247;
t371 = pkin(3) + qJ(5);
t128 = t229 * t371 + t297;
t404 = t251 * t235 + t342 * t236;
t418 = -t225 * pkin(4) + t404;
t136 = t252 * t418;
t63 = -t128 * t250 + t136;
t64 = t252 * t128 + t250 * t418;
t290 = t250 * t64 + t252 * t63;
t51 = -pkin(5) * t225 + t136 + (-pkin(8) * t229 - t128) * t250;
t332 = t229 * t252;
t53 = pkin(8) * t332 + t64;
t37 = -t253 * t53 + t376 * t51;
t38 = t253 * t51 + t376 * t53;
t419 = t278 * t37 + t279 * t38;
t417 = t279 * t270 + t278 * t408;
t415 = mrSges(7,2) / 0.2e1;
t249 = t252 ^ 2;
t414 = t249 / 0.2e1;
t387 = t270 / 0.2e1;
t248 = t250 ^ 2;
t320 = t248 + t249;
t413 = m(6) * t320;
t372 = mrSges(4,3) + mrSges(5,1);
t412 = -Ifges(4,4) - Ifges(5,6);
t411 = mrSges(6,3) * t320;
t187 = mrSges(7,1) * t279 + mrSges(7,2) * t278;
t293 = t250 * mrSges(6,1) + t252 * mrSges(6,2);
t405 = t187 + t293;
t347 = t252 * mrSges(6,1);
t350 = t250 * mrSges(6,2);
t403 = -t350 / 0.2e1 + t347 / 0.2e1;
t314 = t278 ^ 2 + t279 ^ 2;
t179 = -mrSges(6,3) * t225 * t250 + t229 * mrSges(6,1);
t336 = t225 * t252;
t181 = -t229 * mrSges(6,2) + mrSges(6,3) * t336;
t127 = t225 * t371 + t272;
t137 = t252 * t149;
t61 = -t127 * t250 + t137;
t62 = t252 * t127 + t250 * t149;
t402 = -m(6) * (t250 * t61 - t252 * t62) - t250 * t179 + t252 * t181;
t105 = -mrSges(7,2) * t229 + t157 * mrSges(7,3);
t107 = mrSges(7,1) * t229 - t160 * mrSges(7,3);
t380 = -t279 / 0.2e1;
t280 = -t278 * t107 / 0.2e1 + t105 * t380;
t352 = t279 * mrSges(7,3);
t306 = t352 / 0.2e1;
t354 = t278 * mrSges(7,3);
t307 = -t354 / 0.2e1;
t263 = t157 * t306 + t160 * t307 + t280;
t12 = t263 - t407;
t401 = t12 * qJD(1);
t334 = t279 * t160;
t335 = t278 * t157;
t381 = t278 / 0.2e1;
t265 = (-t335 / 0.2e1 - t334 / 0.2e1) * mrSges(7,3) + t105 * t381 + t107 * t380;
t281 = mrSges(7,1) * t387 + t408 * t415;
t14 = t265 - t281;
t400 = t14 * qJD(1);
t398 = 0.2e1 * t225;
t397 = 0.2e1 * t229;
t396 = m(5) / 0.2e1;
t395 = m(5) / 0.4e1;
t393 = m(6) / 0.2e1;
t390 = m(4) * pkin(2);
t388 = t408 / 0.2e1;
t301 = t342 * pkin(2);
t244 = -t301 - pkin(3);
t237 = -qJ(5) + t244;
t370 = -pkin(8) + t237;
t213 = t370 * t250;
t214 = t370 * t252;
t164 = t213 * t376 + t253 * t214;
t386 = -t164 / 0.2e1;
t364 = Ifges(7,4) * t278;
t189 = -Ifges(7,2) * t279 + t364;
t385 = t189 / 0.2e1;
t222 = Ifges(7,4) * t279;
t191 = Ifges(7,1) * t278 - t222;
t384 = t191 / 0.2e1;
t383 = -t225 / 0.2e1;
t382 = t225 / 0.2e1;
t379 = -t250 / 0.2e1;
t378 = t252 / 0.2e1;
t375 = m(5) * t229;
t374 = pkin(2) * t251;
t367 = Ifges(6,4) * t250;
t366 = Ifges(6,4) * t252;
t365 = Ifges(7,4) * t160;
t363 = Ifges(6,5) * t250;
t362 = Ifges(6,6) * t252;
t217 = t229 * mrSges(5,2);
t218 = t229 * mrSges(4,1);
t163 = -t253 * t213 + t214 * t376;
t238 = qJ(4) + t374;
t234 = pkin(5) * t250 + t238;
t305 = -t352 / 0.2e1;
t317 = t390 / 0.2e1;
t337 = t225 * t238;
t257 = (t229 * t244 - t337) * t396 + (t229 * t237 * t320 - t337) * t393 + (t163 * t408 + t164 * t270 - t225 * t234) * t391 + (-t225 * t251 - t229 * t342) * t317 + t408 * t307 + t270 * t305 + t405 * t383 - t229 * t411 / 0.2e1;
t106 = mrSges(7,2) * t225 + mrSges(7,3) * t408;
t108 = -mrSges(7,1) * t225 - mrSges(7,3) * t270;
t178 = pkin(3) * t229 + t297;
t333 = t229 * t250;
t180 = -t225 * mrSges(6,1) - mrSges(6,3) * t333;
t182 = t225 * mrSges(6,2) + mrSges(6,3) * t332;
t260 = t178 * t396 + (-t250 * t63 + t252 * t64) * t393 + (t278 * t38 - t279 * t37) * t391 + t106 * t381 + t108 * t380 + t180 * t379 + t182 * t378 + t254 * t317;
t6 = t225 * t421 - t217 + t218 - t257 + t260;
t361 = qJD(1) * t6;
t104 = t225 * t302 + t404;
t125 = -Ifges(6,6) * t225 + (t252 * Ifges(6,2) + t367) * t229;
t349 = t250 * Ifges(6,1);
t126 = -Ifges(6,5) * t225 + (t349 + t366) * t229;
t294 = -t347 + t350;
t175 = t294 * t225;
t176 = t294 * t229;
t282 = Ifges(7,5) * t387 + Ifges(7,6) * t388;
t50 = pkin(5) * t229 + t137 + (-pkin(8) * t225 - t127) * t250;
t52 = pkin(8) * t336 + t62;
t35 = -t253 * t52 + t376 * t50;
t36 = t253 * t50 + t376 * t52;
t69 = Ifges(7,2) * t157 + t229 * Ifges(7,6) + t365;
t70 = Ifges(7,4) * t270 + Ifges(7,2) * t408 - Ifges(7,6) * t225;
t155 = Ifges(7,4) * t157;
t71 = Ifges(7,1) * t160 + t229 * Ifges(7,5) + t155;
t72 = Ifges(7,1) * t270 + Ifges(7,4) * t408 - Ifges(7,5) * t225;
t357 = t160 * mrSges(7,2);
t359 = t157 * mrSges(7,1);
t77 = t357 - t359;
t78 = t356 - t358;
t1 = -t177 * t217 + t245 * t218 + (mrSges(4,2) * t247 + (t362 + t363 + t412) * t229 + (Ifges(6,2) * t414 + Ifges(5,3) - Ifges(7,3) - Ifges(4,1) - Ifges(5,2) + Ifges(4,2) - Ifges(6,3) + (t366 + t349 / 0.2e1) * t250) * t225 + t282) * t229 + (-Ifges(7,5) * t160 / 0.2e1 - Ifges(7,6) * t157 / 0.2e1 + t250 * t126 / 0.2e1 + t125 * t378 + t177 * mrSges(5,3) - t245 * mrSges(4,2) + (-t363 / 0.2e1 - t362 / 0.2e1 - t412) * t225) * t225 + m(6) * (-t149 * t418 + t61 * t63 + t62 * t64) + t418 * t176 + t71 * t387 + t69 * t388 - t149 * t175 + m(7) * (t103 * t104 + t35 * t37 + t36 * t38) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t254 + (m(4) * t245 + mrSges(4,1) * t225) * pkin(2)) * t254 + ((Ifges(3,1) - Ifges(3,2)) * t254 - pkin(1) * mrSges(3,2) + Ifges(3,4) * t377) * t377 + t423 * t178 + t103 * t77 + t104 * t78 + t38 * t105 + t36 * t106 + t37 * t107 + t35 * t108 + t157 * t70 / 0.2e1 + t160 * t72 / 0.2e1 + t63 * t179 + t61 * t180 + t64 * t181 + t62 * t182;
t360 = t1 * qJD(1);
t355 = t225 * mrSges(5,1);
t325 = Ifges(7,5) * t157 - Ifges(7,6) * t160;
t76 = t160 * mrSges(7,1) + mrSges(7,2) * t157;
t79 = -Ifges(7,2) * t160 + t155;
t80 = Ifges(7,1) * t157 - t365;
t4 = t104 * t76 + t229 * t325 / 0.2e1 - t36 * t107 + t35 * t105 + (-t69 / 0.2e1 - t36 * mrSges(7,3) + t80 / 0.2e1) * t160 + (t79 / 0.2e1 - t35 * mrSges(7,3) + t71 / 0.2e1) * t157;
t345 = t4 * qJD(1);
t292 = t250 * t62 + t252 * t61;
t328 = t252 * t179;
t331 = t250 * t181;
t5 = t270 * t105 + t408 * t107 + (t229 * t372 + t328 + t331) * t229 + (t225 * t372 - t175 - t77) * t225 + m(7) * (-t104 * t225 + t270 * t36 + t35 * t408) + m(6) * (-t225 * t418 + t229 * t292) + (m(5) + m(4)) * (t198 * t229 - t225 * t404);
t344 = t5 * qJD(1);
t271 = t417 * t391;
t276 = m(7) * t417;
t25 = t271 + t375 + t276 / 0.2e1 + (t248 / 0.2e1 + t414) * m(6) * t397;
t341 = qJD(1) * t25;
t16 = m(7) * (t270 * t35 - t36 * t408) + t270 * t107 - t408 * t105 + (-t402 - t423) * t229;
t338 = t16 * qJD(1);
t330 = t250 * t182;
t327 = t252 * t180;
t321 = -Ifges(7,5) * t279 - Ifges(7,6) * t278;
t185 = -mrSges(7,1) * t278 + mrSges(7,2) * t279;
t318 = t185 * qJD(6);
t316 = t394 + t392;
t315 = -m(7) / 0.4e1 - m(6) / 0.4e1;
t308 = t354 / 0.2e1;
t300 = t333 / 0.2e1;
t295 = -t413 / 0.2e1;
t84 = mrSges(5,3) + m(7) * t234 + 0.4e1 * (m(6) / 0.4e1 + t395) * t238 + t405;
t258 = -m(5) * t404 / 0.2e1 + t418 * t394 + (t163 * t270 - t164 * t408 + t104) * t392 + t359 / 0.2e1 - t357 / 0.2e1 + t270 * t308 - t408 * t306 + t403 * t225;
t262 = t404 * t396 + t290 * t393 + t419 * t391 + t108 * t381 + t279 * t106 / 0.2e1 + t330 / 0.2e1 + t327 / 0.2e1;
t9 = t258 + t262;
t289 = -qJD(1) * t9 + qJD(2) * t84;
t285 = qJD(1) * t76 - qJD(2) * t185;
t277 = m(7) * (t334 + t335);
t18 = m(7) * (t157 * t36 - t160 * t35) - t160 * t107 + t157 * t105 + t402 * t225;
t275 = t18 * qJD(1) + qJD(4) * t313;
t269 = -t314 * t391 + t295;
t86 = t269 + t316;
t274 = qJD(1) * t313 + t86 * qJD(2);
t48 = -t277 / 0.2e1 + (-t413 / 0.4e1 + t315) * t398;
t273 = t48 * qJD(1);
t188 = -Ifges(7,2) * t278 - t222;
t190 = -Ifges(7,1) * t279 - t364;
t24 = -t234 * t185 - (t188 / 0.2e1 + t384) * t279 - (t385 - t190 / 0.2e1) * t278;
t259 = -(t69 / 0.4e1 - t80 / 0.4e1) * t278 - (t79 / 0.4e1 + t71 / 0.4e1) * t279 + (-t163 * mrSges(7,3) / 0.2e1 + t188 / 0.4e1 + t191 / 0.4e1) * t157 + (mrSges(7,3) * t386 - t189 / 0.4e1 + t190 / 0.4e1) * t160 - t104 * t185 / 0.2e1 + t163 * t105 / 0.2e1 + t107 * t386 + t229 * t321 / 0.4e1 + t234 * t76 / 0.2e1;
t264 = Ifges(7,3) * t382 - t37 * mrSges(7,1) / 0.2e1 + t38 * t415 - t282;
t3 = t259 + t264;
t268 = -t3 * qJD(1) - t24 * qJD(2);
t261 = t157 * t305 + t160 * t308 + t292 * t394 + (t164 * t157 - t163 * t160 - t278 * t35 - t279 * t36) * t391 - t331 / 0.2e1 - t328 / 0.2e1 + t280;
t11 = t229 * t403 + t261 - t424;
t46 = m(7) * (-t163 * t278 - t164 * t279) + mrSges(7,3) * t314 - t237 * t413 + t411;
t266 = t11 * qJD(1) + t46 * qJD(2);
t85 = t269 - t316;
t54 = qJD(5) * t313;
t49 = t277 / 0.2e1 + t413 * t382 + t315 * t398;
t26 = -t276 / 0.2e1 + t229 * t295 - t375 / 0.2e1 + t271 + (t395 + t413 / 0.4e1) * t397;
t15 = t265 + t281;
t13 = t263 + t407;
t10 = mrSges(6,2) * t300 - mrSges(6,1) * t332 / 0.2e1 + t261 + t424;
t8 = -t258 + t262 - t355;
t7 = t257 + t260;
t2 = t259 - t264;
t17 = [qJD(2) * t1 + qJD(3) * t5 + qJD(4) * t16 + qJD(5) * t18 + qJD(6) * t4, t7 * qJD(3) + t8 * qJD(4) + t10 * qJD(5) + t2 * qJD(6) + t360 + (t70 * t380 + t72 * t381 + t126 * t378 + t125 * t379 + (Ifges(6,5) * t252 + Ifges(7,5) * t278 - Ifges(6,6) * t250 - Ifges(7,6) * t279) * t383 + (-m(6) * t238 - t293) * t149 + t408 * t385 + (m(5) * t244 - t342 * t390 - mrSges(4,1) + mrSges(5,2)) * t404 + m(7) * (t103 * t234 + t163 * t37 + t164 * t38) + (mrSges(3,2) * pkin(7) - Ifges(3,6)) * t254 + Ifges(3,5) * t377 + (-Ifges(6,2) * t250 + t366) * t332 / 0.2e1 + (Ifges(6,1) * t252 - t367) * t300 - t244 * t355 - mrSges(3,1) * t312 + (m(6) * t290 + t327 + t330) * t237 + (-t238 * mrSges(5,1) - mrSges(4,3) * t374 + Ifges(5,5) - Ifges(4,6)) * t229 + (mrSges(4,3) * t301 + Ifges(5,4) - Ifges(4,5)) * t225 + t270 * t384 - t419 * mrSges(7,3) - t290 * mrSges(6,3) + (-m(5) * t238 - t251 * t390 - t421) * t198 + t163 * t108 + t164 * t106 + t103 * t187 + t234 * t78 + t238 * t176) * qJD(2), t7 * qJD(2) + qJD(3) * t425 + t26 * qJD(4) + t49 * qJD(5) + t13 * qJD(6) + t344, t8 * qJD(2) + t26 * qJD(3) + qJD(4) * t425 + t15 * qJD(6) + t338 + t54, t10 * qJD(2) + t49 * qJD(3) + t275, t345 + t2 * qJD(2) + t13 * qJD(3) + t15 * qJD(4) + (-mrSges(7,1) * t36 - mrSges(7,2) * t35 + t325) * qJD(6); -qJD(3) * t6 - qJD(4) * t9 + qJD(5) * t11 + qJD(6) * t3 - t360, qJD(4) * t84 + qJD(5) * t46 + qJD(6) * t24, -t361, qJD(5) * t85 + t289, t85 * qJD(4) + t266 (-mrSges(7,1) * t164 - mrSges(7,2) * t163 + t321) * qJD(6) - t268; qJD(2) * t6 - qJD(4) * t25 - qJD(5) * t48 + qJD(6) * t12 - t344, t361, 0, -t341, -t273, t318 + t401; qJD(2) * t9 + qJD(3) * t25 + qJD(6) * t14 - t338 + t54, qJD(5) * t86 - t289, t341, 0, t274, -t187 * qJD(6) + t400; -t11 * qJD(2) + t48 * qJD(3) + t76 * qJD(6) - t275, -t86 * qJD(4) - t266 - t318, t273, -t274, 0, t285; -qJD(2) * t3 - qJD(3) * t12 - qJD(4) * t14 - qJD(5) * t76 - t345, t185 * qJD(5) + t268, -t401, -t400, -t285, 0;];
Cq  = t17;

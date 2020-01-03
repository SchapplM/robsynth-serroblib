% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:11
% EndTime: 2019-12-31 21:19:18
% DurationCPUTime: 3.40s
% Computational Cost: add. (7735->344), mult. (15239->457), div. (0->0), fcn. (15324->6), ass. (0->227)
t222 = sin(qJ(2));
t201 = (-pkin(7) - pkin(6)) * t222;
t356 = cos(qJ(2));
t294 = t356 * pkin(6);
t202 = t356 * pkin(7) + t294;
t221 = sin(qJ(3));
t355 = cos(qJ(3));
t147 = -t355 * t201 + t202 * t221;
t383 = t221 * t201 + t355 * t202;
t411 = (t147 * t221 + t355 * t383) * pkin(2);
t220 = sin(qJ(5));
t190 = t221 * t222 - t355 * t356;
t394 = -t190 * pkin(4) + t383;
t410 = t220 * t394;
t223 = cos(qJ(5));
t409 = t223 * t394;
t210 = -t356 * pkin(2) - pkin(1);
t191 = t221 * t356 + t355 * t222;
t314 = t191 * qJ(4);
t241 = t210 - t314;
t118 = t190 * pkin(3) + t241;
t384 = m(5) * t118 - mrSges(5,2) * t190 - mrSges(5,3) * t191;
t264 = Ifges(6,5) * t223 - Ifges(6,6) * t220;
t246 = t190 * t264;
t357 = t223 / 0.2e1;
t331 = t220 * Ifges(6,2);
t348 = Ifges(6,4) * t223;
t197 = -t331 + t348;
t303 = t223 * t197;
t349 = Ifges(6,4) * t220;
t198 = t223 * Ifges(6,1) - t349;
t307 = t220 * t198;
t395 = t303 / 0.2e1 + t307 / 0.2e1;
t265 = Ifges(6,2) * t223 + t349;
t73 = -Ifges(6,6) * t190 + t265 * t191;
t266 = Ifges(6,1) * t220 + t348;
t75 = -Ifges(6,5) * t190 + t266 * t191;
t236 = -t246 / 0.2e1 - t220 * t73 / 0.2e1 + t75 * t357 + (Ifges(5,4) - Ifges(4,5)) * t190 + (Ifges(5,5) - Ifges(4,6) + t395) * t191;
t214 = t220 * mrSges(6,1);
t215 = t223 * mrSges(6,2);
t382 = t215 + t214;
t99 = pkin(4) * t191 + t147;
t397 = t99 * t382;
t398 = t383 * mrSges(5,2);
t399 = t383 * mrSges(4,1);
t401 = t147 * mrSges(5,3);
t402 = t147 * mrSges(4,2);
t408 = t236 + t398 + t402 - t397 - t399 - t401;
t407 = t397 / 0.2e1 - t398 / 0.2e1 + t399 / 0.2e1 + t401 / 0.2e1 - t402 / 0.2e1;
t359 = t220 / 0.2e1;
t406 = -t223 / 0.2e1;
t404 = qJ(4) * t99;
t403 = t394 * t99;
t352 = pkin(2) * t221;
t207 = qJ(4) + t352;
t400 = t207 * t99;
t218 = t220 ^ 2;
t219 = t223 ^ 2;
t300 = t218 + t219;
t396 = t300 * mrSges(6,3) + mrSges(4,1);
t289 = mrSges(5,3) + t382;
t346 = Ifges(6,6) * t223;
t347 = Ifges(6,5) * t220;
t393 = Ifges(4,4) + Ifges(5,6) - t347 / 0.2e1 - t346 / 0.2e1;
t392 = -pkin(3) * t383 - qJ(4) * t147;
t293 = t355 * pkin(2);
t209 = -t293 - pkin(3);
t391 = -t207 * t147 + t209 * t383;
t328 = t223 * mrSges(6,1);
t332 = t220 * mrSges(6,2);
t195 = t328 - t332;
t115 = t195 * t190;
t116 = t195 * t191;
t313 = t191 * t220;
t121 = -t190 * mrSges(6,1) - mrSges(6,3) * t313;
t312 = t191 * t223;
t123 = t190 * mrSges(6,2) + mrSges(6,3) * t312;
t366 = pkin(3) + pkin(8);
t76 = t366 * t190 + t241;
t39 = -t220 * t76 + t223 * t99;
t40 = t220 * t99 + t223 * t76;
t389 = t99 * t115 + (t393 * t190 + t73 * t357 + t75 * t359) * t190 - t394 * t116 + t118 * (-mrSges(5,2) * t191 + mrSges(5,3) * t190) + t39 * t121 + t40 * t123 + t210 * (mrSges(4,1) * t191 - mrSges(4,2) * t190);
t128 = pkin(3) * t191 + qJ(4) * t190;
t187 = t191 * pkin(8);
t87 = t128 + t187;
t46 = t223 * t87 + t410;
t323 = t46 * t220;
t45 = -t220 * t87 + t409;
t324 = t45 * t223;
t260 = t323 + t324;
t373 = m(6) / 0.2e1;
t386 = t260 * t373;
t305 = t223 * t121;
t308 = t220 * t123;
t385 = t308 / 0.2e1 + t305 / 0.2e1;
t381 = t305 + t308;
t380 = m(6) / 0.4e1 + m(5) / 0.4e1;
t379 = Ifges(6,6) * t312 / 0.2e1 + Ifges(6,5) * t313 / 0.2e1 - Ifges(6,3) * t190 / 0.2e1;
t378 = t265 * t359 + t266 * t406;
t175 = -(-m(6) - m(5)) * qJ(4) + t289;
t377 = 0.2e1 * m(6);
t376 = 2 * qJD(3);
t375 = m(5) / 0.2e1;
t371 = -mrSges(6,1) / 0.2e1;
t370 = mrSges(6,2) / 0.2e1;
t369 = Ifges(6,3) / 0.2e1;
t333 = t191 * Ifges(6,6);
t72 = t265 * t190 + t333;
t368 = -t72 / 0.4e1;
t365 = t394 / 0.2e1;
t364 = t394 / 0.4e1;
t117 = t382 * t190;
t363 = t117 / 0.2e1;
t316 = t190 * t220;
t120 = t191 * mrSges(6,1) - mrSges(6,3) * t316;
t362 = -t120 / 0.2e1;
t315 = t190 * t223;
t122 = -t191 * mrSges(6,2) + mrSges(6,3) * t315;
t361 = t122 / 0.2e1;
t173 = Ifges(6,4) * t315;
t360 = -t173 / 0.4e1;
t353 = m(5) * t383;
t217 = t222 * pkin(2);
t350 = mrSges(6,3) * t190;
t119 = t128 + t217;
t326 = t223 * t72;
t334 = t191 * Ifges(6,5);
t74 = Ifges(6,1) * t316 + t173 + t334;
t329 = t220 * t74;
t227 = t326 / 0.2e1 + t329 / 0.2e1 - t393 * t191 + (-Ifges(4,1) + Ifges(5,3) + Ifges(4,2) - Ifges(5,2) - Ifges(6,3)) * t190;
t77 = t119 + t187;
t41 = -t220 * t77 + t409;
t42 = t223 * t77 + t410;
t1 = t41 * t120 + t42 * t122 + m(6) * (t39 * t41 + t40 * t42 - t403) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t222 + (m(4) * t210 + mrSges(4,1) * t190) * pkin(2)) * t222 + (mrSges(4,2) * t217 + t227) * t191 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2)) * t222 + Ifges(3,4) * t356) * t356 + t389 + t384 * t119;
t345 = t1 * qJD(1);
t336 = t190 * mrSges(5,1);
t335 = t191 * mrSges(5,1);
t330 = t220 * t42;
t327 = t223 * t41;
t3 = t45 * t120 + t46 * t122 + m(6) * (t39 * t45 + t40 * t46 - t403) + t227 * t191 + t384 * t128 + t389;
t325 = t3 * qJD(1);
t184 = -t307 / 0.2e1;
t262 = t220 * t40 + t223 * t39;
t5 = -t39 * t122 + t40 * t120 - t394 * t117 + (t72 * t359 - t191 * t264 / 0.2e1 + (t331 * t357 + t184) * t190 + t262 * mrSges(6,3) + (t173 + t74) * t406) * t190;
t322 = t5 * qJD(1);
t321 = qJ(4) * t116;
t304 = t223 * t122;
t309 = t220 * t120;
t12 = (m(6) * (t220 * t39 - t223 * t40) - t304 + t309 - t384) * t191;
t320 = qJD(1) * t12;
t240 = (t215 / 0.2e1 + t214 / 0.2e1) * t191;
t247 = t304 / 0.2e1 - t309 / 0.2e1;
t274 = t218 / 0.2e1 + t219 / 0.2e1;
t16 = t274 * t350 + t240 - t247;
t317 = t16 * qJD(1);
t311 = t207 * t116;
t310 = t207 * t195;
t298 = mrSges(6,3) * t327;
t297 = mrSges(6,3) * t330;
t296 = mrSges(6,3) * t324;
t295 = mrSges(6,3) * t323;
t292 = t207 * t335;
t290 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t288 = t355 * mrSges(4,2);
t286 = -t336 / 0.2e1;
t285 = t335 / 0.2e1;
t284 = t360 - t74 / 0.4e1;
t283 = mrSges(5,2) * t352 + t289 * t293;
t282 = t355 * t207;
t281 = t195 * t365;
t280 = -t308 / 0.2e1;
t278 = -t305 / 0.2e1;
t276 = t366 * t362;
t275 = t366 * t361;
t272 = t300 * t221;
t271 = t190 * t293;
t270 = t184 - t303 / 0.2e1 + t378;
t269 = mrSges(6,3) * t274;
t194 = qJ(4) * t195;
t268 = -t194 / 0.2e1 - t310 / 0.2e1;
t263 = -t346 - t347;
t261 = t327 + t330;
t205 = -pkin(8) + t209;
t28 = -t283 + (t288 - m(6) * (t205 * t272 + t282) - m(5) * (t209 * t221 + t282)) * pkin(2) + t396 * t352;
t225 = (t391 + t411) * t375 + (-t400 + (t262 * t221 + t355 * t394) * pkin(2)) * t373 - t311 / 0.2e1 - t292 / 0.2e1 + t209 * t286 - t296 / 0.2e1 - t295 / 0.2e1 - t115 * t293 / 0.2e1 + t285 * t352 - mrSges(5,1) * t271 / 0.2e1 + (t223 * t120 + t220 * t122) * t352 / 0.2e1 + (t385 + t386) * t205 - t407;
t226 = -m(5) * t392 / 0.2e1 - m(6) * (-t261 * t366 - t404) / 0.2e1 + t321 / 0.2e1 + pkin(3) * t286 + qJ(4) * t285 - t366 * t280 - t366 * t278 + t298 / 0.2e1 + t297 / 0.2e1 + t407;
t4 = t225 + t226;
t259 = t4 * qJD(1) - t28 * qJD(2);
t48 = t270 + t310;
t242 = t198 / 0.4e1 + t290 * t223;
t231 = (-0.3e1 / 0.4e1 * t349 + t242) * t190 + t368 - 0.3e1 / 0.4e1 * t333;
t243 = -t197 / 0.4e1 - t290 * t220;
t232 = t243 * t190 - 0.3e1 / 0.4e1 * t334 + t284;
t248 = t207 * t363 + t281;
t253 = t42 * t370 + t41 * t371;
t255 = t205 * t269;
t7 = (t369 - t255) * t190 + (t205 * t362 + t232) * t220 + (t205 * t361 + t231) * t223 + t248 + t253;
t258 = t7 * qJD(1) + t48 * qJD(2);
t124 = 0.4e1 * t380 * t207 + t289;
t250 = t332 / 0.2e1 - t328 / 0.2e1;
t233 = t250 * t190 + t278 + t280;
t13 = (t364 - t330 / 0.4e1 - t327 / 0.4e1) * t377 + t233;
t256 = -qJD(1) * t13 - qJD(2) * t124;
t254 = t366 * t269;
t252 = t46 * t370 + t45 * t371;
t249 = qJ(4) * t363 + t281;
t228 = -t378 + t395;
t237 = t250 * t352;
t31 = -t237 + t228 + t268;
t55 = -t194 + t228;
t9 = (t369 + t254) * t190 + (-t276 + t232) * t220 + (-t275 + t231) * t223 + t249 + t252;
t239 = t9 * qJD(1) - t31 * qJD(2) - t55 * qJD(3);
t18 = (t364 - t323 / 0.4e1 - t324 / 0.4e1) * t377 + t233;
t66 = (-0.1e1 / 0.2e1 + t274) * m(6) * t352 - t175;
t238 = qJD(1) * t18 - qJD(2) * t66 + qJD(3) * t175;
t234 = (-mrSges(5,1) + t250) * t190 + m(6) * t365 + t385;
t229 = t242 * t223 + (-0.3e1 / 0.4e1 * t348 + t243) * t220;
t211 = qJ(4) * t293;
t67 = (t300 * t373 + t375) * t352 + t289 + t380 * (0.4e1 * qJ(4) + 0.2e1 * t352);
t32 = -t237 - t268 + t270;
t17 = t240 + t247 - t300 * t350 / 0.2e1;
t15 = t234 + t353 + t386;
t11 = t261 * t373 + t383 * t375 + t353 / 0.2e1 + t234;
t8 = (-t275 - t333 / 0.4e1 + t368) * t223 + (-t276 - t334 / 0.4e1 + t284) * t220 + (t254 + t229) * t190 + t249 - t252 + t379;
t6 = t191 * t263 / 0.4e1 - t329 / 0.4e1 - t326 / 0.4e1 + t220 * t360 + t247 * t205 + (-t255 + t229) * t190 + t248 - t253 + t379;
t2 = t225 - t226 + t236;
t10 = [qJD(2) * t1 + qJD(3) * t3 + qJD(4) * t12 - qJD(5) * t5, t345 + (-t311 + (pkin(6) * mrSges(3,2) - Ifges(3,6)) * t222 + (-t191 * t352 + t271) * mrSges(4,3) + Ifges(3,5) * t356 - t209 * t336 - mrSges(3,1) * t294 - m(4) * t411 + (m(6) * t261 + t381) * t205 + t408 - m(6) * t400 + m(5) * t391 - t292 - t297 - t298) * qJD(2) + t2 * qJD(3) + t11 * qJD(4) + t6 * qJD(5), t325 + t2 * qJD(2) + (-mrSges(5,1) * t314 + pkin(3) * t336 - t366 * t381 - t295 - t296 - t321 + t408) * qJD(3) + t15 * qJD(4) + t8 * qJD(5) + ((-t260 * t366 - t404) * t373 + t392 * t375) * t376, qJD(2) * t11 + qJD(3) * t15 + qJD(5) * t17 + t320, -t322 + t6 * qJD(2) + t8 * qJD(3) + t17 * qJD(4) + (-mrSges(6,1) * t40 - mrSges(6,2) * t39 + t246) * qJD(5); qJD(3) * t4 + qJD(4) * t13 + qJD(5) * t7 - t345, -qJD(3) * t28 + qJD(4) * t124 + qJD(5) * t48, t67 * qJD(4) + t32 * qJD(5) + ((-pkin(2) * t272 * t366 + t211) * t373 + (-pkin(3) * t352 + t211) * t375) * t376 + t259 + (t283 + (-t396 * t221 - t288) * pkin(2)) * qJD(3), qJD(3) * t67 - t256, t32 * qJD(3) + (-t205 * t382 + t263) * qJD(5) + t258; -qJD(2) * t4 + qJD(4) * t18 + qJD(5) * t9 - t325, -qJD(4) * t66 - qJD(5) * t31 - t259, qJD(4) * t175 - qJD(5) * t55, t238, ((mrSges(6,2) * t366 - Ifges(6,6)) * t223 + (mrSges(6,1) * t366 - Ifges(6,5)) * t220) * qJD(5) + t239; -qJD(2) * t13 - qJD(3) * t18 - qJD(5) * t16 - t320, qJD(3) * t66 + t256, -t238, 0, -qJD(5) * t382 - t317; -qJD(2) * t7 - qJD(3) * t9 + qJD(4) * t16 + t322, qJD(3) * t31 - t258, -t239, t317, 0;];
Cq = t10;

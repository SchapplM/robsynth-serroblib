% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:13
% EndTime: 2019-12-05 16:54:23
% DurationCPUTime: 3.96s
% Computational Cost: add. (4245->409), mult. (11004->570), div. (0->0), fcn. (9999->8), ass. (0->210)
t273 = sin(pkin(5));
t274 = sin(qJ(4));
t276 = sin(qJ(2));
t277 = cos(qJ(4));
t278 = cos(qJ(3));
t279 = cos(qJ(2));
t341 = t278 * t279;
t154 = (-t274 * t341 + t276 * t277) * t273;
t155 = (t274 * t276 + t277 * t341) * t273;
t370 = -qJ(5) - pkin(8);
t234 = t370 * t274;
t237 = t370 * t277;
t259 = -pkin(4) * t277 - pkin(3);
t275 = sin(qJ(3));
t346 = t273 * t279;
t324 = t275 * t346;
t351 = t155 * t277;
t352 = t154 * t274;
t418 = mrSges(6,3) + mrSges(5,3);
t429 = t418 * (t352 / 0.2e1 - t351 / 0.2e1) - m(5) * (-pkin(3) * t324 + (t351 - t352) * pkin(8)) / 0.2e1 - m(6) * (t234 * t154 - t237 * t155 + t259 * t324) / 0.2e1;
t407 = m(5) / 0.2e1;
t406 = m(6) / 0.2e1;
t368 = mrSges(5,1) * t277;
t305 = -mrSges(5,2) * t274 + t368;
t428 = t305 / 0.2e1;
t427 = Ifges(5,1) + Ifges(6,1);
t372 = -Ifges(5,2) - Ifges(6,2);
t426 = t372 * t274;
t345 = t274 * t275;
t335 = mrSges(6,3) * t345;
t357 = t278 * mrSges(6,2);
t219 = -t335 + t357;
t334 = mrSges(5,3) * t345;
t220 = mrSges(5,2) * t278 - t334;
t425 = t220 + t219;
t371 = Ifges(5,6) + Ifges(6,6);
t313 = t371 * t274;
t373 = Ifges(5,5) + Ifges(6,5);
t424 = t373 * t277 - Ifges(4,4) - t313;
t422 = m(5) + m(6);
t421 = 0.2e1 * t406;
t420 = 0.2e1 * t407;
t270 = t275 ^ 2;
t419 = pkin(7) * t270;
t417 = Ifges(5,3) + Ifges(6,3);
t405 = m(6) * pkin(4);
t416 = -mrSges(6,1) - t405;
t415 = t371 * t277;
t263 = t274 * mrSges(6,1);
t264 = t277 * mrSges(6,2);
t414 = t264 + t263;
t246 = pkin(3) * t275 - pkin(8) * t278;
t160 = pkin(7) * t345 + t277 * t246;
t343 = t275 * t277;
t161 = -pkin(7) * t343 + t274 * t246;
t411 = -t160 * t274 + t161 * t277;
t358 = t278 * mrSges(4,2);
t239 = t275 * mrSges(4,1) + t358;
t410 = t264 / 0.2e1 + t263 / 0.2e1;
t344 = t274 * t278;
t100 = -qJ(5) * t344 + t161;
t232 = -pkin(3) * t278 - pkin(8) * t275 - pkin(2);
t216 = t277 * t232;
t336 = pkin(7) * t344;
t156 = t216 - t336;
t342 = t277 * t278;
t157 = pkin(7) * t342 + t232 * t274;
t347 = t273 * t276;
t355 = cos(pkin(5));
t191 = t275 * t347 - t355 * t278;
t192 = t355 * t275 + t278 * t347;
t238 = t274 * mrSges(5,1) + t277 * mrSges(5,2);
t197 = t238 * t275;
t198 = t414 * t278;
t199 = t238 * t278;
t221 = -mrSges(6,2) * t275 - mrSges(6,3) * t344;
t222 = -mrSges(5,2) * t275 - mrSges(5,3) * t344;
t226 = mrSges(5,1) * t275 - mrSges(5,3) * t342;
t378 = pkin(4) * t274;
t323 = pkin(7) + t378;
t227 = t323 * t275;
t228 = t323 * t278;
t394 = -t220 / 0.2e1;
t395 = -t219 / 0.2e1;
t316 = t394 + t395;
t376 = pkin(7) * t278;
t377 = pkin(7) * t275;
t225 = mrSges(6,1) * t275 - mrSges(6,3) * t342;
t389 = t225 / 0.2e1;
t224 = -mrSges(5,1) * t278 - mrSges(5,3) * t343;
t390 = t224 / 0.2e1;
t359 = t278 * mrSges(6,1);
t223 = -mrSges(6,3) * t343 - t359;
t392 = t223 / 0.2e1;
t196 = t414 * t275;
t396 = t196 / 0.2e1;
t306 = -qJ(5) * t343 + t216;
t90 = (-pkin(7) * t274 - pkin(4)) * t278 + t306;
t91 = pkin(4) * t275 - qJ(5) * t342 + t160;
t93 = -t192 * t274 - t277 * t346;
t94 = t192 * t277 - t274 * t346;
t96 = -qJ(5) * t345 + t157;
t408 = (t396 + t197 / 0.2e1) * t192 + (t221 / 0.2e1 + t222 / 0.2e1) * t94 + (t389 + t226 / 0.2e1) * t93 + (t160 * t93 + t161 * t94 + t192 * t377) * t407 + (t100 * t94 + t192 * t227 + t91 * t93) * t406 + ((t274 * t90 - t277 * t96 + t228) * t406 + (t156 * t274 - t157 * t277 + t376) * t407 + t198 / 0.2e1 + t199 / 0.2e1 + t316 * t277 + (t390 + t392) * t274) * t191;
t271 = t277 ^ 2;
t272 = t278 ^ 2;
t404 = -mrSges(5,1) / 0.2e1;
t403 = mrSges(6,1) / 0.2e1;
t402 = mrSges(5,2) / 0.2e1;
t401 = mrSges(6,2) / 0.2e1;
t400 = t93 / 0.2e1;
t95 = t306 - t336;
t399 = t95 / 0.2e1;
t398 = pkin(4) * t91;
t397 = t191 / 0.2e1;
t393 = -t223 / 0.2e1;
t391 = -t224 / 0.2e1;
t388 = t234 / 0.2e1;
t387 = t237 / 0.2e1;
t386 = t238 / 0.2e1;
t385 = -t274 / 0.2e1;
t384 = t274 / 0.2e1;
t382 = t277 / 0.2e1;
t381 = m(6) * t227;
t380 = m(6) * t259;
t379 = m(6) * t275;
t374 = Ifges(6,4) + Ifges(5,4);
t369 = -t90 + t95;
t367 = mrSges(6,3) * t277;
t366 = Ifges(5,4) * t274;
t268 = Ifges(5,4) * t277;
t365 = Ifges(6,4) * t274;
t267 = Ifges(6,4) * t277;
t364 = t274 * t93;
t363 = t274 * t94;
t361 = t277 * t93;
t360 = t277 * t94;
t356 = -t305 - mrSges(4,1);
t301 = t360 - t364;
t11 = (t192 - t301) * t191 * t422;
t354 = t11 * qJD(1);
t348 = t191 * t275;
t12 = (t93 * t154 + t94 * t155 + t191 * t324) * t422 + m(4) * (t192 * t278 - t347 + t348) * t346;
t353 = t12 * qJD(1);
t269 = t274 ^ 2;
t339 = t271 + t269;
t337 = pkin(4) * t343;
t333 = t228 * t406;
t332 = t378 / 0.2e1;
t331 = pkin(8) * t394;
t330 = t402 + t401;
t329 = -mrSges(6,3) / 0.2e1 - mrSges(5,3) / 0.2e1;
t328 = Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1;
t327 = -0.5e1 / 0.4e1 * Ifges(5,4) - 0.5e1 / 0.4e1 * Ifges(6,4);
t326 = Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1;
t325 = m(6) * t369;
t319 = t346 / 0.2e1;
t265 = Ifges(6,5) * t277;
t266 = Ifges(5,5) * t277;
t315 = -t266 / 0.4e1 - t265 / 0.4e1;
t314 = t277 * t374;
t310 = t404 - t405 / 0.2e1;
t309 = t275 * t319;
t308 = t330 * t277;
t256 = mrSges(6,1) * t343;
t307 = -mrSges(6,2) * t345 + t256;
t235 = -mrSges(6,1) * t277 + mrSges(6,2) * t274;
t1 = (t358 / 0.2e1 - t239 / 0.2e1 + (mrSges(4,1) / 0.2e1 + t428 - t235 / 0.2e1) * t275) * t346 + t408 + t429;
t5 = pkin(2) * t239 - t227 * t198 - t228 * t196 - t199 * t377 - t197 * t376 - t90 * t225 - t91 * t223 - t96 * t221 - t100 * t219 - t156 * t226 - t157 * t222 - t160 * t224 - t161 * t220 - m(6) * (t100 * t96 + t227 * t228 + t90 * t91) - m(5) * (t156 * t160 + t157 * t161) + t424 * t272 + (-t424 * t275 + (-Ifges(4,1) + Ifges(4,2) - m(5) * pkin(7) ^ 2 - t427 * t271 + (0.2e1 * t314 + t426) * t274 + t417) * t278) * t275;
t303 = t1 * qJD(1) - t5 * qJD(2);
t298 = t403 - t310;
t285 = t298 * t154 - t330 * t155;
t296 = t325 / 0.2e1 + t393;
t6 = -t191 * t256 / 0.2e1 + (t330 * t274 + t310 * t277) * t348 + (-t329 * t343 - t296 + t390) * t94 + (t329 * t345 + t316) * t93 + t285;
t8 = -t305 * t419 - t196 * t337 - t227 * t307 - t95 * t219 + t96 * t223 - t90 * t335 + t157 * t224 - t96 * t325 + ((-t373 * t278 - t374 * t345) * t274 + (-pkin(4) * t381 + t96 * mrSges(6,3) + t157 * mrSges(5,3) + t275 * t314 - t371 * t278 + (t372 + t427) * t345) * t277) * t275 + (-t220 - t334) * t156;
t302 = -t6 * qJD(1) - t8 * qJD(2);
t28 = (m(6) * (-t274 * t96 - t90 * t277) - t274 * t219 - t277 * t223) * t275;
t34 = (t319 + t363 / 0.2e1 + t361 / 0.2e1) * t379;
t300 = -qJD(1) * t34 + qJD(2) * t28;
t299 = t234 * t274 + t237 * t277;
t141 = -m(6) * t337 - t307;
t200 = -m(6) * t378 - t414;
t297 = qJD(2) * t141 + qJD(3) * t200;
t240 = t277 * Ifges(6,2) + t365;
t241 = t277 * Ifges(5,2) + t366;
t242 = t274 * Ifges(6,1) + t267;
t243 = t274 * Ifges(5,1) + t268;
t16 = pkin(3) * t238 - t259 * t414 - t235 * t378 + (-t242 / 0.2e1 - t267 / 0.2e1 - t243 / 0.2e1 - t268 / 0.2e1) * t277 + (-pkin(4) * t380 + t240 / 0.2e1 + t241 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t274 + (t326 - t328) * t277) * t274;
t282 = mrSges(6,3) * t387 + pkin(3) * t404 - t241 / 0.4e1 - t240 / 0.4e1 + (-Ifges(5,2) / 0.4e1 - Ifges(6,2) / 0.4e1 + t328) * t277 + (t235 / 0.2e1 + t380 / 0.2e1) * pkin(4);
t283 = pkin(3) * t402 + mrSges(6,3) * t388 - t259 * mrSges(6,2) / 0.2e1 - t243 / 0.4e1 - t242 / 0.4e1 - t268 / 0.4e1 - t267 / 0.4e1 + (-Ifges(5,1) / 0.4e1 - Ifges(6,1) / 0.4e1 + t326) * t274;
t287 = t100 * t401 + t160 * t404 + t161 * t402 - t91 * mrSges(6,1) / 0.2e1;
t288 = t227 * t414 / 0.2e1 + t219 * t388 + t259 * t256 / 0.2e1;
t289 = pkin(7) * t386 + (-t271 / 0.2e1 - t269 / 0.2e1) * pkin(8) * mrSges(5,3);
t292 = pkin(8) * t391 + (t399 - t90 / 0.2e1) * mrSges(6,3);
t4 = t223 * t387 + t274 * t331 + (-t225 / 0.2e1 + t196 * t384) * pkin(4) + (t227 * t332 + t90 * t387 - t237 * t399 - t398 / 0.2e1) * m(6) + t292 * t277 + ((-0.3e1 / 0.4e1 * Ifges(5,5) - 0.3e1 / 0.4e1 * Ifges(6,5)) * t277 + t313 + t315) * t278 + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t282 * t277 + (t327 * t277 + t283) * t274 + t289) * t275 + t287 + t288;
t9 = (-t238 / 0.2e1 + t308 + (t403 + mrSges(5,1) / 0.2e1) * t274 - t410) * t191;
t293 = -t9 * qJD(1) + t4 * qJD(2) - t16 * qJD(3);
t286 = m(6) * ((-t234 * t275 + t96) * t277 + (t237 * t275 - t90) * t274);
t21 = (t357 / 0.2e1 + t395) * t277 + (t359 / 0.2e1 + t392) * t274 + t333 - t286 / 0.2e1;
t32 = 0.2e1 * (t192 / 0.4e1 + t364 / 0.4e1 - t360 / 0.4e1) * m(6);
t57 = -m(6) * t299 + t339 * mrSges(6,3);
t291 = -qJD(1) * t32 - qJD(2) * t21 + qJD(3) * t57;
t233 = t346 * t419;
t35 = (-t361 - t363) * t379 / 0.2e1 + m(6) * t309;
t33 = (t192 + t301) * t406;
t22 = t286 / 0.2e1 + t219 * t382 + t223 * t385 + t333 + t410 * t278;
t10 = t414 * t397 + (m(6) * t332 + t298 * t274 + t308 + t386) * t191;
t7 = t348 * t428 + t307 * t397 + t191 * t337 * t406 + t285 + (t369 * t406 + t391 + t393) * t94 + t425 * t400 + t418 * (-t94 * t343 / 0.2e1 + t345 * t400);
t3 = -t287 + ((-Ifges(5,5) / 0.4e1 - Ifges(6,5) / 0.4e1) * t278 + t292) * t277 + (t331 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t278 + (t381 / 0.2e1 + t396) * pkin(4)) * t274 - t296 * t237 + t315 * t278 + pkin(4) * t389 + t398 * t406 + t288 - t371 * t344 / 0.2e1 + t373 * t342 / 0.2e1 + ((t327 * t274 + t282) * t277 + t283 * t274 + t289 + t417 / 0.2e1) * t275;
t2 = (t235 - t305) * t309 - t239 * t346 + t408 - t429;
t13 = [qJD(2) * t12 + qJD(3) * t11, t2 * qJD(3) + t7 * qJD(4) + t35 * qJD(5) + t353 + (((-t278 * mrSges(4,1) + t275 * mrSges(4,2) - mrSges(3,1)) * t276 + (-mrSges(3,2) + (t196 + t197) * t275 + (t270 + t272) * mrSges(4,3)) * t279) * t273 + t233 * t420 + t227 * t324 * t421 + m(4) * (t233 + (pkin(7) * t272 * t279 - pkin(2) * t276) * t273) + (t157 * t420 + t96 * t421 + t425) * t155 + (t156 * t420 + t90 * t421 + t223 + t224) * t154) * qJD(2), t2 * qJD(2) + t10 * qJD(4) + t33 * qJD(5) + t354 + ((t235 + t356) * t192 + (-t418 * t339 + mrSges(4,2)) * t191 + (-t339 * t191 * pkin(8) - pkin(3) * t192) * t420 + (t299 * t191 + t192 * t259) * t421) * qJD(3), t7 * qJD(2) + t10 * qJD(3) + ((-mrSges(5,2) - mrSges(6,2)) * t93 + (-mrSges(5,1) + t416) * t94) * qJD(4), qJD(2) * t35 + qJD(3) * t33; qJD(3) * t1 - qJD(4) * t6 - qJD(5) * t34 - t353, -qJD(3) * t5 - qJD(4) * t8 + qJD(5) * t28, t3 * qJD(4) + t22 * qJD(5) + t303 + (t259 * t198 + t234 * t225 + t228 * t235 - t237 * t221 - pkin(3) * t199 + m(6) * (-t100 * t237 + t228 * t259 + t234 * t91) - t91 * t274 * mrSges(6,3) + t100 * t367 + (m(5) * t411 + t277 * t222 - t274 * t226) * pkin(8) + (pkin(7) * mrSges(4,2) + t373 * t274 - Ifges(4,6) + t415) * t275 + t411 * mrSges(5,3) + (Ifges(4,5) + (-m(5) * pkin(3) + t356) * pkin(7) + (t240 + t241) * t385 + (t427 * t277 - t365 - t366) * t384 + (t242 + t243 + t267 + t268 + t426) * t382) * t278) * qJD(3), t3 * qJD(3) + t302 + (-mrSges(5,1) * t157 - mrSges(5,2) * t156 - mrSges(6,2) * t95 + (-t415 + (mrSges(6,3) * pkin(4) - t373) * t274) * t275 + t416 * t96) * qJD(4), qJD(3) * t22 + t300; -qJD(2) * t1 - qJD(4) * t9 - qJD(5) * t32 - t354, qJD(4) * t4 - qJD(5) * t21 - t303, -qJD(4) * t16 + qJD(5) * t57, t293 + (-mrSges(6,2) * t234 - pkin(4) * t367 - pkin(8) * t368 + t265 + t266 + (mrSges(5,2) * pkin(8) - t371) * t274 - t416 * t237) * qJD(4), t291; qJD(2) * t6 + qJD(3) * t9, -qJD(3) * t4 + qJD(5) * t141 - t302, qJD(5) * t200 - t293, 0, t297; qJD(2) * t34 + qJD(3) * t32, qJD(3) * t21 - qJD(4) * t141 - t300, -qJD(4) * t200 - t291, -t297, 0;];
Cq = t13;

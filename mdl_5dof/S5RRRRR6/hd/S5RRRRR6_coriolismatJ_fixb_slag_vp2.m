% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:14:51
% EndTime: 2020-01-03 12:15:02
% DurationCPUTime: 6.35s
% Computational Cost: add. (17888->277), mult. (34216->359), div. (0->0), fcn. (35872->8), ass. (0->189)
t522 = qJD(3) + qJD(4);
t250 = sin(qJ(4));
t251 = sin(qJ(3));
t254 = cos(qJ(4));
t255 = cos(qJ(3));
t222 = -t250 * t251 + t254 * t255;
t223 = -t250 * t255 - t251 * t254;
t249 = sin(qJ(5));
t253 = cos(qJ(5));
t175 = t222 * t249 - t223 * t253;
t300 = t222 * t253 + t223 * t249;
t487 = Ifges(6,5) * t300 - Ifges(6,6) * t175;
t252 = sin(qJ(2));
t241 = pkin(1) * t252 + pkin(7);
t375 = pkin(8) + t241;
t219 = t375 * t251;
t220 = t375 * t255;
t166 = -t219 * t250 + t220 * t254;
t388 = pkin(9) * t222;
t128 = t166 + t388;
t218 = t223 * pkin(9);
t430 = -t219 * t254 - t220 * t250;
t443 = t430 + t218;
t481 = -t128 * t249 + t253 * t443;
t490 = t481 * mrSges(6,2);
t74 = t128 * t253 + t249 * t443;
t500 = t74 * mrSges(6,1);
t507 = -t500 / 0.2e1 - t490 / 0.2e1;
t515 = t487 + 0.2e1 * t507;
t521 = t515 * qJD(5);
t399 = -pkin(8) - pkin(7);
t235 = t399 * t251;
t236 = t399 * t255;
t185 = t235 * t250 - t236 * t254;
t147 = t185 + t388;
t429 = t235 * t254 + t236 * t250;
t444 = t429 + t218;
t480 = -t147 * t249 + t253 * t444;
t491 = t480 * mrSges(6,2);
t98 = t147 * t253 + t249 * t444;
t499 = t98 * mrSges(6,1);
t508 = -t499 / 0.2e1 - t491 / 0.2e1;
t516 = t487 + 0.2e1 * t508;
t520 = t516 * qJD(5);
t176 = -mrSges(5,1) * t223 + mrSges(5,2) * t222;
t244 = -pkin(3) * t255 - pkin(2);
t256 = cos(qJ(2));
t393 = pkin(1) * t256;
t229 = t244 - t393;
t189 = -pkin(4) * t222 + t244;
t186 = t189 - t393;
t441 = mrSges(6,1) * t175 + mrSges(6,2) * t300;
t462 = t186 * t441;
t518 = t229 * t176 + t462;
t461 = t189 * t441;
t517 = t244 * t176 + t461;
t449 = t429 * mrSges(5,2);
t465 = t185 * mrSges(5,1);
t483 = -t491 - t499;
t514 = -t449 - t465 + t483;
t448 = t430 * mrSges(5,2);
t466 = t166 * mrSges(5,1);
t482 = -t490 - t500;
t513 = -t448 - t466 + t482;
t419 = t507 - t448 / 0.2e1 - t466 / 0.2e1;
t420 = t508 - t449 / 0.2e1 - t465 / 0.2e1;
t451 = Ifges(6,4) * t300 ^ 2 + (-Ifges(6,4) * t175 + (Ifges(6,1) - Ifges(6,2)) * t300) * t175;
t263 = t451 + (t186 / 0.2e1 + t189 / 0.2e1) * t441;
t510 = (t244 / 0.2e1 + t229 / 0.2e1) * t176 + t263;
t354 = t300 * mrSges(6,3);
t495 = t249 * t481 - t253 * t74;
t494 = t249 * t480 - t253 * t98;
t415 = m(5) / 0.2e1;
t390 = pkin(3) * t254;
t242 = pkin(4) + t390;
t321 = t249 * t250;
t205 = -pkin(3) * t321 + t242 * t253;
t209 = (t253 * t254 - t321) * pkin(3);
t486 = -t205 + t209;
t192 = t223 * t393;
t193 = t222 * t393;
t130 = t192 * t253 - t193 * t249;
t131 = t192 * t249 + t193 * t253;
t320 = t250 * t253;
t206 = pkin(3) * t320 + t242 * t249;
t414 = -m(6) / 0.2e1;
t318 = t130 * mrSges(6,1) / 0.2e1 - t131 * mrSges(6,2) / 0.2e1;
t423 = t318 - t193 * mrSges(5,2) / 0.2e1 + t192 * mrSges(5,1) / 0.2e1;
t485 = t423 - (t130 * t205 + t131 * t206) * t414 + (t192 * t254 + t193 * t250) * pkin(3) * t415;
t484 = (Ifges(5,4) * t222 + (-Ifges(5,1) + Ifges(5,2)) * t223) * t222 - Ifges(5,4) * t223 ^ 2;
t389 = pkin(4) * t223;
t477 = m(6) * t389;
t284 = (-t249 * t254 - t320) * pkin(3);
t474 = t206 + t284;
t104 = -mrSges(6,1) * t300 + mrSges(6,2) * t175;
t391 = pkin(3) * t251;
t196 = -t389 + t391;
t471 = t196 * t104 + t484;
t454 = t251 ^ 2 + t255 ^ 2;
t246 = Ifges(4,4) * t255;
t177 = -mrSges(5,1) * t222 - mrSges(5,2) * t223;
t278 = -Ifges(4,4) * t251 + pkin(3) * t177 + (Ifges(4,1) - Ifges(4,2)) * t255;
t428 = t255 * t246 + t278 * t251;
t452 = t451 + t428 + t471;
t308 = -t354 / 0.2e1;
t355 = t175 * mrSges(6,3);
t446 = -(t474 * t480 + t486 * t98) * t414 + t420;
t445 = -(t474 * t481 + t486 * t74) * t414 + t419;
t427 = t454 * t256;
t232 = mrSges(4,1) * t251 + mrSges(4,2) * t255;
t426 = -mrSges(4,1) * t255 + mrSges(4,2) * t251;
t422 = -t104 * t389 + t484;
t119 = t206 * t355;
t417 = t119 / 0.2e1 + t284 * t355 / 0.2e1 + t209 * t308;
t416 = (-t130 * t175 + t131 * t300) * mrSges(6,3) + (t192 * t223 + t193 * t222) * mrSges(5,3);
t413 = m(6) / 0.2e1;
t412 = -pkin(4) / 0.2e1;
t157 = t249 * pkin(4) * t355;
t398 = -t157 / 0.2e1;
t395 = m(5) * t252;
t394 = m(6) * t252;
t392 = pkin(3) * t250;
t363 = pkin(3) * qJD(3);
t362 = pkin(4) * qJD(4);
t347 = t209 * mrSges(6,2);
t346 = t249 * mrSges(6,1);
t136 = t186 * t196;
t217 = t229 * t391;
t243 = -pkin(2) - t393;
t5 = m(5) * t217 + m(6) * t136 + t243 * t232 + t452 + t518;
t345 = t5 * qJD(1);
t261 = t451 + t422;
t6 = -t186 * t477 + t261 + t518;
t344 = t6 * qJD(1);
t197 = t206 * mrSges(6,1);
t61 = mrSges(6,2) * t205 + t197;
t341 = qJD(5) * t61;
t15 = t451 + t462;
t334 = t15 * qJD(1);
t269 = (mrSges(4,3) * t454 - mrSges(3,2)) * t256 + (-mrSges(3,1) + t104 + t177 + t426) * t252;
t21 = m(6) * (t130 * t481 + t131 * t74) + m(5) * (t166 * t193 + t192 * t430) + (t186 * t394 + t229 * t395 + m(4) * (t241 * t427 + t243 * t252) + t269) * pkin(1) + t416;
t328 = t21 * qJD(1);
t316 = t205 * t354;
t314 = t253 * t354;
t313 = t254 * t222 * mrSges(5,3);
t311 = t393 / 0.2e1;
t301 = (-t186 - t189) * t223;
t299 = t253 * t308;
t298 = Ifges(5,5) * t222 + Ifges(5,6) * t223 + t487;
t296 = -t205 / 0.2e1 + t253 * t412;
t148 = t189 * t196;
t230 = t244 * t391;
t257 = (t243 / 0.2e1 - pkin(2) / 0.2e1) * t232 + (t217 + t230) * t415 + (t136 + t148) * t413 + t471 + t510;
t2 = t257 + (mrSges(4,1) * t311 + t278) * t251 + (mrSges(4,2) * t311 + t246) * t255 - t485;
t7 = m(5) * t230 + m(6) * t148 - pkin(2) * t232 + t452 + t517;
t294 = t2 * qJD(1) + t7 * qJD(2);
t293 = -t157 + t298;
t10 = -t189 * t477 + t261 + t517;
t258 = t422 + t510;
t286 = m(6) * (t130 * t253 + t131 * t249);
t4 = t258 + 0.2e1 * (-t286 / 0.4e1 + m(6) * t301 / 0.4e1) * pkin(4) - t423;
t292 = t4 * qJD(1) + t10 * qJD(2);
t18 = t451 + t461;
t9 = t263 - t318;
t291 = qJD(1) * t9 + qJD(2) * t18;
t289 = -t197 / 0.2e1 + t346 * t412;
t283 = t495 * t413;
t282 = t494 * t413;
t262 = t398 + t316 / 0.2e1 + pkin(4) * t299 + t417;
t12 = pkin(4) * t283 + t262 + t419 - t445;
t14 = pkin(4) * t282 + t262 + t420 - t446;
t207 = mrSges(6,1) * t284;
t276 = -mrSges(5,1) * t392 - mrSges(5,2) * t390 + t207 - t347;
t60 = -m(6) * (t205 * t284 + t206 * t209) - t276;
t281 = -qJD(1) * t12 - qJD(2) * t14 - qJD(3) * t60;
t280 = t61 * qJD(3);
t228 = (mrSges(6,2) * t253 + t346) * pkin(4);
t49 = -t207 / 0.2e1 + (t209 / 0.2e1 + t296) * mrSges(6,2) + t289;
t274 = -qJD(3) * t49 + qJD(4) * t228;
t271 = t205 * t308 + t298 + t398 - t417;
t270 = mrSges(5,3) * t223 * t392 + Ifges(4,5) * t255 - Ifges(4,6) * t251 - t119 + t298 - t316;
t221 = t228 * qJD(5);
t50 = -t347 / 0.2e1 + t207 / 0.2e1 + t296 * mrSges(6,2) + t289;
t13 = (t282 + t299) * pkin(4) + t271 + t420 + t446;
t11 = (t283 + t299) * pkin(4) + t271 + t419 + t445;
t8 = t263 + t318;
t3 = t258 + t423 + (t301 * t413 + t286 / 0.2e1) * pkin(4);
t1 = t257 - t232 * t393 / 0.2e1 + t428 + t485;
t16 = [qJD(2) * t21 + qJD(3) * t5 + qJD(4) * t6 + qJD(5) * t15, t1 * qJD(3) + t3 * qJD(4) + t8 * qJD(5) + t328 + (0.2e1 * (t130 * t480 + t131 * t98) * t413 + 0.2e1 * (t185 * t193 + t192 * t429) * t415 + (t189 * t394 + t244 * t395 + m(4) * (-pkin(2) * t252 + pkin(7) * t427) + t269) * pkin(1) + t416) * qJD(2), t345 + t1 * qJD(2) + (m(6) * (-t205 * t74 + t206 * t481) + t270 + t426 * t241 + t513) * qJD(3) + t11 * qJD(4) + t521 + (-t313 + m(5) * (-t166 * t254 + t250 * t430)) * t363, t344 + t3 * qJD(2) + t11 * qJD(3) + (t293 + t513) * qJD(4) + t521 + (m(6) * t495 - t314) * t362, t334 + t8 * qJD(2) + (t487 + t482) * qJD(5) + t522 * t515; qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t9 - t328, qJD(3) * t7 + qJD(4) * t10 + qJD(5) * t18, (m(6) * (-t205 * t98 + t206 * t480) + t270 + t426 * pkin(7) + t514) * qJD(3) + t13 * qJD(4) + t520 + (-t313 + m(5) * (-t185 * t254 + t250 * t429)) * t363 + t294, t13 * qJD(3) + (t293 + t514) * qJD(4) + t520 + (m(6) * t494 - t314) * t362 + t292, (t487 + t483) * qJD(5) + t291 + t522 * t516; -qJD(2) * t2 - qJD(4) * t12 - t345, -qJD(4) * t14 - t294, -qJD(4) * t60 - t341, (m(6) * (t209 * t249 + t253 * t284) * pkin(4) + t276) * qJD(4) + t50 * qJD(5) + t281, qJD(4) * t50 - t280 - t341; -qJD(2) * t4 + qJD(3) * t12 - t344, qJD(3) * t14 - t292, qJD(5) * t49 - t281, -t221, -t221 - t274; -qJD(2) * t9 - t334, -t291, -qJD(4) * t49 + t280, t274, 0;];
Cq = t16;

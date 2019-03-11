% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:06
% EndTime: 2019-03-09 03:18:15
% DurationCPUTime: 5.28s
% Computational Cost: add. (10034->453), mult. (19556->589), div. (0->0), fcn. (20754->6), ass. (0->232)
t249 = sin(pkin(9));
t250 = cos(pkin(9));
t252 = sin(qJ(3));
t371 = cos(qJ(3));
t226 = t249 * t252 - t250 * t371;
t251 = sin(qJ(5));
t247 = t251 ^ 2;
t253 = cos(qJ(5));
t248 = t253 ^ 2;
t295 = t248 / 0.2e1 + t247 / 0.2e1;
t267 = (mrSges(6,3) + mrSges(7,3)) * t295;
t411 = t267 * t226;
t228 = t249 * t371 + t252 * t250;
t241 = pkin(5) * t251 + qJ(4);
t382 = pkin(3) + pkin(8);
t397 = -qJ(6) - t382;
t231 = t397 * t251;
t232 = t397 * t253;
t274 = t231 * t251 + t232 * t253;
t324 = t247 + t248;
t291 = t324 * t228;
t342 = qJ(4) * t226;
t145 = pkin(3) * t228 + t342;
t370 = m(5) * t145;
t319 = t370 / 0.2e1;
t391 = -m(7) / 0.2e1;
t393 = -m(6) / 0.2e1;
t410 = -t267 * t228 - t319 - (-t291 * t382 - t342) * t393 - (-t226 * t241 + t228 * t274) * t391;
t373 = -t251 / 0.2e1;
t372 = -t253 / 0.2e1;
t403 = Ifges(6,5) + Ifges(7,5);
t401 = Ifges(6,6) + Ifges(7,6);
t275 = t231 * t253 - t232 * t251;
t366 = t228 * pkin(5);
t308 = -pkin(2) * t250 - pkin(1);
t269 = -qJ(4) * t228 + t308;
t68 = t226 * t382 + t269;
t292 = qJ(6) * t226 + t68;
t361 = pkin(7) + qJ(2);
t233 = t361 * t249;
t234 = t361 * t250;
t152 = t371 * t233 + t234 * t252;
t87 = pkin(4) * t228 + t152;
t82 = t253 * t87;
t41 = -t251 * t292 + t82;
t33 = t41 + t366;
t348 = t251 * t87;
t42 = t253 * t292 + t348;
t281 = t251 * t42 + t253 * t33;
t339 = t226 * t251;
t137 = t228 * mrSges(7,1) - mrSges(7,3) * t339;
t338 = t226 * t253;
t141 = -t228 * mrSges(7,2) + mrSges(7,3) * t338;
t329 = t137 * t372 + t141 * t373;
t407 = (t226 * t275 - t281) * t391 - t329;
t353 = t228 * mrSges(6,1);
t138 = -mrSges(6,3) * t339 + t353;
t380 = t138 / 0.2e1;
t377 = -t226 / 0.2e1;
t406 = t228 / 0.2e1;
t405 = m(7) * (t33 * t251 - t253 * t42);
t365 = mrSges(6,1) + mrSges(7,1);
t364 = mrSges(6,2) + mrSges(7,2);
t404 = mrSges(5,3) - mrSges(4,2);
t402 = -Ifges(5,6) - Ifges(4,4);
t400 = Ifges(6,3) + Ifges(7,3);
t389 = m(7) * pkin(5);
t399 = -mrSges(7,1) - t389;
t317 = Ifges(7,4) / 0.2e1 + Ifges(6,4) / 0.2e1;
t398 = t317 * t253;
t244 = t253 * mrSges(7,1);
t293 = -t251 * mrSges(7,2) + t244;
t315 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t358 = Ifges(7,4) * t251;
t282 = Ifges(7,2) * t253 + t358;
t69 = t228 * Ifges(7,6) + t226 * t282;
t360 = Ifges(6,4) * t251;
t283 = Ifges(6,2) * t253 + t360;
t71 = t228 * Ifges(6,6) + t226 * t283;
t394 = t315 * t228 + t69 / 0.2e1 + t71 / 0.2e1;
t392 = m(6) / 0.2e1;
t390 = m(7) / 0.2e1;
t388 = mrSges(6,2) / 0.2e1;
t387 = -mrSges(7,2) / 0.2e1;
t386 = mrSges(7,2) / 0.2e1;
t385 = -mrSges(7,3) / 0.2e1;
t153 = -t252 * t233 + t234 * t371;
t89 = -t226 * pkin(4) + t153;
t83 = t253 * t89;
t84 = t228 * t382 + t342;
t36 = -pkin(5) * t226 + t83 + (-qJ(6) * t228 - t84) * t251;
t384 = -t36 / 0.2e1;
t383 = t89 / 0.2e1;
t367 = pkin(5) * t253;
t307 = -pkin(4) - t367;
t58 = t226 * t307 + t153;
t381 = m(7) * t58;
t337 = t228 * t251;
t139 = -t226 * mrSges(7,1) - mrSges(7,3) * t337;
t379 = -t139 / 0.2e1;
t352 = t228 * mrSges(6,2);
t142 = mrSges(6,3) * t338 - t352;
t378 = t142 / 0.2e1;
t235 = t251 * mrSges(7,1) + t253 * mrSges(7,2);
t376 = t235 / 0.2e1;
t369 = m(7) * t226;
t368 = m(7) * t241;
t363 = mrSges(4,3) + mrSges(5,1);
t52 = t251 * t89 + t253 * t84;
t359 = Ifges(6,4) * t253;
t357 = Ifges(7,4) * t253;
t121 = t293 * t226;
t347 = t253 * mrSges(6,1);
t350 = t251 * mrSges(6,2);
t285 = t347 - t350;
t122 = t285 * t226;
t327 = t141 + t142;
t328 = t137 + t138;
t47 = -t251 * t68 + t82;
t48 = t253 * t68 + t348;
t5 = (t226 * t363 + t121 + t122) * t226 + (t363 * t228 + t327 * t251 + t328 * t253) * t228 + m(7) * (-t226 * t58 + t228 * t281) + m(6) * (-t226 * t89 + (t251 * t48 + t253 * t47) * t228) + (m(3) * qJ(2) + mrSges(3,3)) * (t249 ^ 2 + t250 ^ 2) + (m(5) + m(4)) * (t152 * t228 - t153 * t226);
t356 = qJD(1) * t5;
t209 = t228 * mrSges(5,2);
t210 = t228 * mrSges(4,1);
t140 = -t226 * mrSges(6,1) - mrSges(6,3) * t337;
t336 = t228 * t253;
t143 = t226 * mrSges(7,2) + mrSges(7,3) * t336;
t144 = t226 * mrSges(6,2) + mrSges(6,3) * t336;
t44 = qJ(6) * t336 + t52;
t51 = -t251 * t84 + t83;
t258 = (t379 - t140 / 0.2e1) * t251 + (t144 / 0.2e1 + t143 / 0.2e1) * t253 + t319 + (-t251 * t51 + t253 * t52) * t392 + (-t251 * t36 + t253 * t44) * t390;
t346 = t253 * mrSges(6,2);
t351 = t251 * mrSges(6,1);
t8 = t210 - t209 + (t351 / 0.2e1 + t346 / 0.2e1 + t376 + t404) * t226 + t258 - t410;
t355 = qJD(1) * t8;
t120 = pkin(3) * t226 + t269;
t123 = t293 * t228;
t124 = t285 * t228;
t146 = -mrSges(5,2) * t226 - mrSges(5,3) * t228;
t310 = (t357 + t359 + (Ifges(6,1) + Ifges(7,1)) * t251) * t406 + t403 * t377;
t204 = Ifges(7,4) * t338;
t73 = Ifges(7,1) * t339 + t228 * Ifges(7,5) + t204;
t205 = Ifges(6,4) * t338;
t75 = Ifges(6,1) * t339 + t228 * Ifges(6,5) + t205;
t311 = t73 / 0.2e1 + t75 / 0.2e1;
t312 = (t282 + t283) * t406 + t401 * t377;
t316 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t57 = t228 * t307 - t152;
t1 = t308 * t210 + t145 * t146 + t36 * t137 + t51 * t138 + t33 * t139 + t47 * t140 + t44 * t141 + t52 * t142 + t42 * t143 + t48 * t144 - t57 * t121 + t87 * t122 - t58 * t123 - t89 * t124 + m(6) * (t47 * t51 + t48 * t52 - t87 * t89) + m(7) * (t33 * t36 + t42 * t44 + t57 * t58) + (-t308 * mrSges(4,2) + t310 * t251 + t312 * t253 + (-t251 * t316 - t253 * t315 - t402) * t226) * t226 + (t402 * t228 + t394 * t253 + (t228 * t316 + t311) * t251 + (-Ifges(5,2) + Ifges(5,3) + Ifges(4,2) - Ifges(4,1) - t400) * t226) * t228 + (t226 * mrSges(5,3) - t209 + t370) * t120;
t354 = t1 * qJD(1);
t125 = mrSges(7,1) * t339 + mrSges(7,2) * t338;
t126 = t226 * (t346 + t351);
t127 = -Ifges(7,2) * t339 + t204;
t128 = -Ifges(6,2) * t339 + t205;
t238 = t253 * Ifges(7,1) - t358;
t129 = t226 * t238;
t239 = t253 * Ifges(6,1) - t360;
t130 = t226 * t239;
t202 = Ifges(7,5) * t338;
t203 = Ifges(6,5) * t338;
t309 = m(7) * (-t33 + t41);
t4 = t58 * t125 + t89 * t126 - t48 * t138 + t41 * t141 + t47 * t142 + (t203 / 0.2e1 + t202 / 0.2e1) * t228 + (t309 - t137) * t42 + ((t128 / 0.2e1 + t127 / 0.2e1 - t33 * mrSges(7,3) - t47 * mrSges(6,3) + t311) * t253 + (t130 / 0.2e1 + t129 / 0.2e1 - t42 * mrSges(7,3) - t48 * mrSges(6,3) + (-t121 + t381) * pkin(5) - t394) * t251) * t226;
t345 = t4 * qJD(1);
t288 = t309 / 0.2e1;
t268 = -t137 / 0.2e1 + t288;
t320 = t389 / 0.2e1;
t290 = mrSges(6,1) / 0.2e1 + t320;
t271 = mrSges(7,1) / 0.2e1 + t290;
t298 = -t141 / 0.2e1 - t142 / 0.2e1;
t318 = t388 + t386;
t6 = (t228 * t318 + t298) * t253 + t411 + (t228 * t271 - t268 + t380) * t251;
t344 = t6 * qJD(1);
t314 = -t41 / 0.2e1 + t33 / 0.2e1;
t303 = t336 / 0.2e1;
t326 = mrSges(7,1) * t303 + t337 * t387;
t9 = (-t352 / 0.2e1 - t298) * t251 + (t353 / 0.2e1 + t137 / 0.2e1 + t380 + (t366 / 0.2e1 + t314) * m(7)) * t253 + t326;
t343 = t9 * qJD(1);
t14 = (-m(5) * t120 - t146 - t327 * t253 + t328 * t251 + t405 + m(6) * (t251 * t47 - t253 * t48)) * t228;
t341 = qJD(1) * t14;
t332 = t253 * t141;
t335 = t251 * t137;
t20 = (t332 - t335 - t405) * t226;
t340 = qJD(1) * t20;
t289 = (m(7) / 0.4e1 + m(6) / 0.4e1) * t291;
t294 = m(7) * t324;
t31 = -0.2e1 * t289 + 0.2e1 * (-m(5) / 0.2e1 - m(6) * t324 / 0.4e1 - t294 / 0.4e1) * t228;
t331 = t31 * qJD(1);
t61 = 0.2e1 * (t247 / 0.4e1 + t248 / 0.4e1 + 0.1e1 / 0.4e1) * t369;
t330 = t61 * qJD(1);
t323 = qJD(5) * t251;
t322 = qJD(5) * t253;
t223 = (-0.1e1 / 0.2e1 - t295) * m(7);
t321 = t223 * qJD(3);
t304 = t337 / 0.2e1;
t236 = -t251 * Ifges(7,2) + t357;
t237 = -t251 * Ifges(6,2) + t359;
t297 = t236 / 0.2e1 + t237 / 0.2e1;
t296 = t238 / 0.2e1 + t239 / 0.2e1;
t286 = mrSges(7,2) * t304 + t390 * t57 - mrSges(7,1) * t336 / 0.2e1;
t279 = t251 * t52 + t253 * t51;
t24 = -qJ(4) * t285 - t235 * t367 - t241 * t293 + (-t317 * t251 + t296) * t251 + (-pkin(5) * t368 + t398 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.2e1 + Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1) * t251 + t297) * t253;
t260 = -t237 / 0.4e1 - t236 / 0.4e1 + t231 * t385 + (-Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1) * t251 - t398 + (t368 / 0.2e1 + t376) * pkin(5);
t261 = t268 * t231 + qJ(4) * t126 / 0.2e1 + t232 * t141 / 0.2e1 + t241 * t125 / 0.2e1 + t58 * t244 / 0.2e1;
t262 = t314 * mrSges(7,3) - t127 / 0.4e1 - t128 / 0.4e1 - t73 / 0.4e1 - t75 / 0.4e1 + t382 * t380 + t58 * t387 - t89 * mrSges(6,2) / 0.2e1;
t263 = (-t121 / 0.2e1 + t381 / 0.2e1) * pkin(5) + t129 / 0.4e1 + t130 / 0.4e1 - t69 / 0.4e1 - t71 / 0.4e1 - t382 * t378 + mrSges(6,1) * t383;
t265 = mrSges(7,1) * t384 + t44 * t386 - t51 * mrSges(6,1) / 0.2e1 + t52 * t388;
t266 = t239 / 0.4e1 + t238 / 0.4e1 + t232 * t385 + (-Ifges(6,2) / 0.4e1 - Ifges(7,2) / 0.4e1) * t253;
t270 = t295 * t382 * mrSges(6,3);
t3 = (m(7) * t384 + t379) * pkin(5) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t270) * t226 + ((-0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t228 + t266 * t226 + t263) * t253 + ((-0.3e1 / 0.4e1 * Ifges(6,5) - 0.3e1 / 0.4e1 * Ifges(7,5)) * t228 + t260 * t226 + t262) * t251 + t261 + t265;
t278 = t3 * qJD(1) - t24 * qJD(3);
t18 = t286 + t407;
t86 = -m(7) * t274 + mrSges(7,3) * t324;
t277 = -qJD(1) * t18 + qJD(3) * t86;
t256 = m(6) * t383 + (-t228 * t275 + t58) * t390 + t364 * t339 / 0.2e1 - t365 * t338 / 0.2e1;
t259 = t279 * t393 + (t251 * t44 + t253 * t36) * t391 + (t143 + t144) * t373 + (t139 + t140) * t372;
t13 = t256 + t259;
t135 = t368 + mrSges(5,3) + t364 * t253 + t365 * t251 + (m(6) + m(5)) * qJ(4);
t273 = qJD(1) * t13 + qJD(3) * t135;
t242 = m(7) * t367;
t213 = -t242 - t293;
t66 = -t339 * t389 - t125;
t272 = qJD(1) * t66 + qJD(3) * t213;
t222 = t324 * t391 + t390;
t60 = -t369 / 0.2e1 + t226 * t294 / 0.2e1;
t30 = -0.2e1 * t289 + (m(6) + m(7)) * t291 / 0.2e1;
t19 = t286 - t407;
t12 = m(5) * t153 - t226 * mrSges(5,1) + t256 - t259;
t11 = -t126 / 0.2e1 + t235 * t377 + t258 + t410;
t10 = t253 * t288 + t138 * t372 + t142 * t373 + (-t350 / 0.2e1 + t290 * t253) * t228 + t326 + t329;
t7 = t251 * t288 + t332 / 0.2e1 - t335 / 0.2e1 + t138 * t373 + t253 * t378 + (t251 * t271 + t253 * t318) * t228 - t411;
t2 = (t251 * t260 + t253 * t266 + t270) * t226 + ((-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t228 + t263) * t253 + t36 * t320 + pkin(5) * t139 / 0.2e1 + ((-Ifges(6,5) / 0.4e1 - Ifges(7,5) / 0.4e1) * t228 + t262) * t251 + t261 - t265 + t400 * t377 + t403 * t304 + t401 * t303;
t15 = [qJD(2) * t5 + qJD(3) * t1 + qJD(4) * t14 + qJD(5) * t4 + qJD(6) * t20, qJD(3) * t11 + qJD(4) * t30 + qJD(5) * t10 + qJD(6) * t60 + t356, t11 * qJD(2) + t12 * qJD(4) + t2 * qJD(5) + t19 * qJD(6) + t354 + (-qJ(4) * t124 - t241 * t123 + t232 * t139 + t231 * t143 + t57 * t235 + (-mrSges(6,2) * t87 - t51 * mrSges(6,3) - t36 * mrSges(7,3) - t140 * t382 + t310) * t253 + (-mrSges(6,1) * t87 - t52 * mrSges(6,3) - t44 * mrSges(7,3) - t144 * t382 - t312) * t251 + (pkin(3) * mrSges(5,1) + t251 * t315 - t253 * t316 + Ifges(5,4) - Ifges(4,5)) * t226 + 0.2e1 * (-qJ(4) * t87 - t279 * t382) * t392 + 0.2e1 * (t231 * t44 + t232 * t36 + t241 * t57) * t390 + (-qJ(4) * mrSges(5,1) + t251 * t296 + t253 * t297 + Ifges(5,5) - Ifges(4,6)) * t228 + (-m(5) * pkin(3) - mrSges(4,1) + mrSges(5,2)) * t153 + (-m(5) * qJ(4) - t404) * t152) * qJD(3), qJD(2) * t30 + qJD(3) * t12 + qJD(5) * t7 + t341, t10 * qJD(2) + t2 * qJD(3) + t7 * qJD(4) + t345 + (-mrSges(6,1) * t48 - mrSges(6,2) * t47 - mrSges(7,2) * t41 + t202 + t203 + (-mrSges(7,3) * t367 - t251 * t401) * t226 + t399 * t42) * qJD(5), qJD(2) * t60 + qJD(3) * t19 + t340; qJD(3) * t8 + qJD(4) * t31 - qJD(5) * t9 + qJD(6) * t61 - t356, 0, t355, t331, -t343 + (t251 * t364 - t242 - t244 - t347) * qJD(5), t330; -qJD(2) * t8 + qJD(4) * t13 + qJD(5) * t3 - qJD(6) * t18 - t354, -t355, qJD(4) * t135 - qJD(5) * t24 + qJD(6) * t86, qJD(6) * t222 + t273 (-mrSges(7,2) * t232 + t231 * t399) * qJD(5) + (mrSges(6,2) * t382 - t401) * t322 + (mrSges(6,1) * t382 + mrSges(7,3) * pkin(5) - t403) * t323 + t278, qJD(4) * t222 + t277; -qJD(2) * t31 - qJD(3) * t13 - qJD(5) * t6 - t341, -t331, qJD(6) * t223 - t273, 0, -t344 - t364 * t322 + (-t365 - t389) * t323, t321; qJD(2) * t9 - qJD(3) * t3 + qJD(4) * t6 + qJD(6) * t66 - t345, t343, t213 * qJD(6) - t278, t344, 0, t272; -qJD(2) * t61 + qJD(3) * t18 - qJD(5) * t66 - t340, -t330, -qJD(4) * t223 - qJD(5) * t213 - t277, -t321, -t272, 0;];
Cq  = t15;

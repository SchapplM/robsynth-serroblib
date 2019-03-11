% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP4
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:20
% EndTime: 2019-03-09 03:11:25
% DurationCPUTime: 3.63s
% Computational Cost: add. (5329->433), mult. (10199->566), div. (0->0), fcn. (8021->6), ass. (0->226)
t222 = sin(qJ(5));
t219 = t222 ^ 2;
t224 = cos(qJ(5));
t220 = t224 ^ 2;
t296 = t219 + t220;
t388 = (mrSges(7,2) + mrSges(6,3)) * t296;
t317 = t224 * mrSges(6,2);
t326 = t222 * mrSges(6,1);
t167 = t317 + t326;
t213 = t224 * mrSges(7,3);
t325 = t222 * mrSges(7,1);
t376 = t213 - t325;
t387 = t167 - t376;
t223 = sin(qJ(3));
t355 = t223 / 0.2e1;
t382 = Ifges(7,4) + Ifges(6,5);
t206 = sin(pkin(9)) * pkin(1) + pkin(7);
t340 = pkin(4) + t206;
t225 = cos(qJ(3));
t274 = t220 / 0.2e1 + t219 / 0.2e1;
t270 = t274 * mrSges(7,2);
t375 = t270 + (-t219 / 0.2e1 + t274) * mrSges(6,3);
t385 = t375 * t225;
t384 = -Ifges(6,6) / 0.2e1;
t383 = -t223 / 0.2e1;
t352 = t225 / 0.2e1;
t289 = m(7) * t355;
t381 = Ifges(7,2) + Ifges(6,3);
t380 = -Ifges(5,6) - Ifges(4,4);
t379 = t223 * (m(7) / 0.4e1 + m(6) / 0.4e1);
t300 = t224 * t225;
t155 = -t223 * mrSges(6,2) - mrSges(6,3) * t300;
t320 = t223 * mrSges(7,3);
t156 = -mrSges(7,2) * t300 + t320;
t378 = t155 + t156;
t311 = qJ(6) * t224;
t264 = pkin(5) * t222 - t311;
t377 = -m(7) * t264 - t387;
t215 = Ifges(7,5) * t224;
t174 = -Ifges(7,1) * t222 + t215;
t373 = -Ifges(7,3) * t222 - t215;
t212 = m(7) * qJ(6) + mrSges(7,3);
t110 = Ifges(7,4) * t223 + t225 * t174;
t333 = Ifges(6,4) * t224;
t268 = Ifges(6,1) * t222 + t333;
t112 = Ifges(6,5) * t223 - t268 * t225;
t147 = t340 * t225;
t287 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t226 = -pkin(3) - pkin(8);
t284 = -cos(pkin(9)) * pkin(1) - pkin(2);
t251 = -qJ(4) * t223 + t284;
t104 = t226 * t225 + t251;
t146 = t340 * t223;
t51 = -t104 * t222 + t146 * t224;
t46 = -pkin(5) * t223 - t51;
t312 = qJ(6) * t222;
t347 = pkin(5) * t224;
t169 = t312 + t347;
t63 = t169 * t225 + t147;
t372 = t287 * t223 + t147 * mrSges(6,2) + t46 * mrSges(7,2) + t110 / 0.2e1 + t112 / 0.2e1 - t51 * mrSges(6,3) - t63 * mrSges(7,3);
t371 = 2 * qJD(3);
t370 = -m(6) / 0.2e1;
t369 = m(6) / 0.2e1;
t368 = -m(7) / 0.2e1;
t367 = m(7) / 0.2e1;
t366 = -pkin(5) / 0.2e1;
t365 = -mrSges(7,3) / 0.2e1;
t364 = Ifges(7,6) / 0.2e1;
t301 = t224 * t104;
t307 = t222 * t146;
t52 = t301 + t307;
t363 = t52 / 0.2e1;
t313 = qJ(4) * t225;
t168 = -t223 * pkin(3) + t313;
t152 = pkin(8) * t223 - t168;
t59 = t147 * t224 - t152 * t222;
t56 = -pkin(5) * t225 - t59;
t362 = t56 / 0.2e1;
t361 = t63 / 0.2e1;
t360 = pkin(5) * mrSges(7,2);
t359 = -qJ(4) / 0.2e1;
t358 = t147 / 0.2e1;
t304 = t222 * t225;
t321 = t223 * mrSges(7,1);
t154 = -mrSges(7,2) * t304 - t321;
t357 = t154 / 0.2e1;
t356 = t222 / 0.2e1;
t354 = -t224 / 0.2e1;
t353 = t224 / 0.2e1;
t351 = -t226 / 0.2e1;
t162 = qJ(4) + t264;
t350 = m(7) * t162;
t349 = m(7) * t169;
t346 = t52 * mrSges(6,1);
t345 = t59 * mrSges(6,1);
t344 = mrSges(6,1) + mrSges(7,1);
t343 = mrSges(5,2) - mrSges(4,1);
t342 = mrSges(6,2) - mrSges(7,3);
t341 = -mrSges(5,3) + mrSges(4,2);
t303 = t223 * qJ(6);
t45 = t52 + t303;
t339 = t45 - t52;
t338 = t46 + t51;
t337 = m(7) * qJD(6);
t336 = mrSges(6,3) * t222;
t334 = Ifges(6,4) * t222;
t332 = Ifges(7,5) * t222;
t331 = Ifges(6,6) * t223;
t330 = Ifges(7,6) * t223;
t328 = t162 * mrSges(7,1);
t327 = t162 * mrSges(7,3);
t324 = t222 * mrSges(6,2);
t323 = t222 * mrSges(7,3);
t322 = t223 * mrSges(6,1);
t319 = t224 * mrSges(6,1);
t318 = t224 * mrSges(7,1);
t316 = t225 * mrSges(7,1);
t141 = t264 * t225;
t306 = t222 * t154;
t241 = t306 / 0.2e1 + t378 * t353;
t260 = t51 * t222 - t52 * t224;
t261 = -t46 * t222 - t45 * t224;
t302 = t223 * t224;
t279 = t302 / 0.2e1;
t7 = t141 * t289 + (mrSges(6,2) * t279 + t376 * t383 + (-t260 + t261) * t368 + t241 + t385) * t225;
t315 = t7 * qJD(1);
t228 = (t338 * t222 + t339 * t224) * t368 - t241;
t230 = t264 * t367 + t325 / 0.2e1 + t317 / 0.2e1 - t213 / 0.2e1;
t8 = -t385 + (t230 + t326) * t223 + t228;
t314 = t8 * qJD(1);
t244 = t222 * (mrSges(6,3) * t304 + t322);
t139 = -pkin(3) * t225 + t251;
t272 = -m(5) * t139 - t225 * mrSges(5,2) + mrSges(5,3) * t223;
t12 = (-m(6) * t260 - m(7) * t261 + t378 * t224 - t244 - t272 + t306) * t223;
t310 = qJD(1) * t12;
t163 = t318 + t323;
t140 = t163 * t225;
t21 = m(7) * (t223 * t45 + t63 * t304) + t140 * t304 + t223 * t156;
t309 = qJD(1) * t21;
t308 = t162 * t225;
t305 = t222 * t223;
t299 = t225 * t376;
t298 = t226 * t223;
t60 = t222 * t147 + t224 * t152;
t297 = t296 * t298;
t295 = qJD(3) * t224;
t294 = qJD(3) * t225;
t293 = qJD(5) * t225;
t291 = m(7) * t305;
t290 = t222 * t337;
t288 = mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t286 = t364 + t384;
t285 = t363 - t45 / 0.2e1;
t280 = -t302 / 0.2e1;
t265 = -Ifges(7,3) * t224 + t332;
t266 = Ifges(6,2) * t224 + t334;
t278 = Ifges(7,6) * t352 + t225 * t384 + t265 * t355 + t266 * t383;
t277 = -t155 / 0.2e1 - t156 / 0.2e1;
t173 = -Ifges(6,2) * t222 + t333;
t276 = t373 / 0.2e1 + t173 / 0.2e1;
t175 = Ifges(7,1) * t224 + t332;
t177 = Ifges(6,1) * t224 - t334;
t275 = t175 / 0.2e1 + t177 / 0.2e1;
t273 = m(5) * t206 + mrSges(5,1);
t271 = t360 - t382;
t269 = 0.2e1 * t296 * t379;
t164 = t319 - t324;
t106 = -t265 * t225 + t330;
t108 = -t266 * t225 + t331;
t239 = t45 * mrSges(7,2) - t106 / 0.2e1 + t108 / 0.2e1 - t147 * mrSges(6,1) - t63 * mrSges(7,1);
t253 = -t59 * mrSges(6,3) + (-t174 + t268) * t355 + t382 * t352;
t53 = qJ(6) * t225 + t60;
t62 = (-t169 - t340) * t223;
t3 = t62 * t140 + t56 * t154 + t60 * t155 + t53 * t156 + t272 * t168 + m(7) * (t45 * t53 + t46 * t56 + t62 * t63) + m(6) * (-t146 * t147 + t51 * t59 + t52 * t60) + (t345 - t139 * mrSges(5,2) + t284 * mrSges(4,1) + t380 * t223 + (t52 * mrSges(6,3) - t286 * t223 + t239) * t224 + t372 * t222) * t223 + (-t139 * mrSges(5,3) + t284 * mrSges(4,2) - t52 * mrSges(6,2) + t51 * mrSges(6,1) - t46 * mrSges(7,1) + t45 * mrSges(7,3) + (-t146 * mrSges(6,1) + t278) * t224 + (t146 * mrSges(6,2) - t253) * t222 + (-Ifges(5,3) + Ifges(5,2) + Ifges(4,1) - Ifges(4,2) + t381) * t223 + (-t287 * t222 + t286 * t224 - t380) * t225) * t225;
t246 = t318 / 0.2e1 + t323 / 0.2e1;
t258 = t222 * t60 + t224 * t59;
t259 = t222 * t53 - t224 * t56;
t6 = ((t357 + t321 / 0.2e1) * t224 + (t288 * t223 + t277) * t222 + (t222 * t45 - t224 * t46 + t62) * t368 + (t222 * t52 + t224 * t51 - t146) * t370) * t223 + (t280 * t336 - t140 / 0.2e1 + t246 * t225 + (-t259 + t63) * t368 + (t147 - t258) * t370) * t225;
t263 = t3 * qJD(1) - t6 * qJD(2);
t142 = t373 * t225;
t143 = t225 * t173;
t144 = t225 * t175;
t145 = t225 * t177;
t205 = Ifges(6,6) * t304;
t5 = t52 * t154 + m(7) * (-t141 * t63 + t46 * t52) - t141 * t140 + (t205 / 0.2e1 - t346) * t223 + ((-t330 / 0.2e1 + t144 / 0.2e1 + t145 / 0.2e1 + t239) * t222 + (t142 / 0.2e1 + t143 / 0.2e1 - t372) * t224) * t225 + (m(7) * t45 + t378) * t51;
t262 = t5 * qJD(1) - t7 * qJD(2);
t35 = 0.4e1 * (0.1e1 - t296) * t225 * t379;
t257 = -t6 * qJD(1) + t35 * qJD(2);
t232 = (-t224 * t63 + (t298 + t308) * t222) * t367 + t140 * t354;
t249 = m(7) * t362 - t316 / 0.2e1;
t19 = (mrSges(7,2) * t223 + t299 / 0.2e1) * t222 - t232 + t249;
t64 = (-t376 + t350) * t224;
t256 = qJD(1) * t19 + qJD(3) * t64;
t24 = t320 + 0.2e1 * (t303 / 0.2e1 + t301 / 0.4e1 + t307 / 0.4e1 - t52 / 0.4e1) * m(7);
t254 = qJD(1) * t24 + qJD(5) * t212;
t252 = m(6) * t358 + m(7) * t361;
t250 = t338 * t367 + t357;
t248 = mrSges(7,2) * t305 - t316;
t247 = mrSges(7,2) * t302 + t225 * mrSges(7,3);
t18 = qJ(4) * t164 - t169 * t376 + (t163 + t349) * t162 + (t174 / 0.2e1 - t268 / 0.2e1 - t276) * t224 + (-t265 / 0.2e1 + t266 / 0.2e1 - t275) * t222;
t214 = Ifges(7,6) * t224;
t227 = (-t162 * t141 + t169 * t63) * t367 + t141 * t376 / 0.2e1 + t164 * t358 + t169 * t140 / 0.2e1 + t223 * t214 / 0.4e1 + t163 * t361;
t229 = (-pkin(5) * t56 + qJ(6) * t53) * t368 + t53 * t365 + mrSges(7,1) * t362 - t345 / 0.2e1 + t60 * mrSges(6,2) / 0.2e1;
t231 = (t339 * t367 - t277) * t226 + t106 / 0.4e1 - t108 / 0.4e1 - t144 / 0.4e1 - t145 / 0.4e1;
t236 = mrSges(6,2) * t359 + t327 / 0.2e1 - t177 / 0.4e1 - t175 / 0.4e1 + t266 / 0.4e1 - t265 / 0.4e1;
t237 = mrSges(6,1) * t359 - t328 / 0.2e1 + t268 / 0.4e1 - t174 / 0.4e1 + t173 / 0.4e1 + t373 / 0.4e1;
t238 = t143 / 0.4e1 + t142 / 0.4e1 - t112 / 0.4e1 - t110 / 0.4e1 + (-t51 / 0.2e1 - t46 / 0.2e1) * mrSges(7,2);
t2 = (mrSges(7,1) * t366 - Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1 + qJ(6) * t365 + (t274 * mrSges(6,3) + t270) * t226) * t225 + (t250 * t226 + (t360 / 0.2e1 - 0.3e1 / 0.4e1 * Ifges(6,5) - 0.3e1 / 0.4e1 * Ifges(7,4) + mrSges(6,1) * t351) * t223 + (t336 * t351 + t237) * t225 + t238) * t222 + ((-0.3e1 / 0.4e1 * Ifges(6,6) + t364) * t223 + (-t303 / 0.2e1 + t285) * mrSges(7,2) + t236 * t225 + t231) * t224 + t227 + t229;
t240 = (mrSges(7,1) / 0.2e1 + mrSges(6,1) / 0.2e1) * t224 + t288 * t222;
t29 = (-t164 / 0.2e1 - t163 / 0.2e1 + (t347 / 0.2e1 + t312 / 0.2e1 - t169 / 0.2e1) * m(7) + t240) * t223;
t243 = t2 * qJD(1) - t29 * qJD(2) + t18 * qJD(3);
t233 = t258 * t370 + t259 * t368;
t11 = t233 + t252;
t48 = (t368 + t370) * t223 + t269;
t55 = t317 - t213 + t350 + mrSges(5,3) + t344 * t222 + (m(6) + m(5)) * qJ(4);
t242 = qJD(1) * t11 - qJD(2) * t48 + qJD(3) * t55;
t157 = (m(7) * t226 - mrSges(7,2)) * t222;
t39 = m(5) * t223 + m(6) * t355 + t269 + t289;
t30 = t169 * t289 + (-t324 / 0.2e1 + t319 / 0.2e1 + t349 / 0.2e1 + t246) * t223 + (t164 + t163) * t355;
t23 = (t52 + 0.2e1 * t303) * t367 + m(7) * t363 + t156;
t20 = -t299 * t356 + t232 + t249;
t10 = t248 * t354 + (t225 * mrSges(6,1) - mrSges(6,3) * t305) * t353 + (t240 + t273) * t225 - t233 + t252 + (-t225 * mrSges(6,2) + mrSges(6,3) * t302 + t247) * t356;
t9 = -t244 / 0.2e1 + (t326 / 0.2e1 + t230) * t223 - t228 + t352 * t388;
t4 = -qJD(3) * t6 - qJD(5) * t7;
t1 = (t237 * t222 + t236 * t224 + t375 * t226) * t225 + t248 * t366 + qJ(6) * t247 / 0.2e1 + Ifges(6,6) * t279 + Ifges(7,6) * t280 + ((-Ifges(6,5) / 0.4e1 - Ifges(7,4) / 0.4e1) * t223 + (-t322 / 0.2e1 + t250) * t226 + t238) * t222 + (-t331 / 0.4e1 + t285 * mrSges(7,2) + t231) * t224 + t227 - t229 + t381 * t352 + t382 * t305 / 0.2e1;
t13 = [qJD(3) * t3 - qJD(4) * t12 + qJD(5) * t5 + qJD(6) * t21, t4, t10 * qJD(4) + t1 * qJD(5) + t20 * qJD(6) + (t56 * mrSges(7,2) + t253) * t295 + ((-qJ(4) * t146 + t258 * t226) * t369 + (t162 * t62 + t259 * t226) * t367) * t371 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + Ifges(4,5) + (-m(5) * pkin(3) + t343) * t206 + (t344 * t226 + t287) * t224 + (-t342 * t226 + t286) * t222) * t294 + t263 + (-t146 * t167 - t62 * t376 + (-t53 * mrSges(7,2) - t60 * mrSges(6,3) + t278) * t222 + (Ifges(5,5) - Ifges(4,6) + (t276 - t328) * t224 + (t275 - t327) * t222 + t341 * t206 + (-t164 - t273) * qJ(4)) * t223) * qJD(3), qJD(3) * t10 + qJD(5) * t9 - t310, t1 * qJD(3) + t9 * qJD(4) + (t205 - t346 + (-m(7) * pkin(5) - mrSges(7,1)) * t52 + (-mrSges(6,2) + t212) * t51) * qJD(5) + t23 * qJD(6) + ((qJ(6) * mrSges(7,2) - Ifges(7,6)) * t222 + t271 * t224) * t293 + t262, qJD(3) * t20 + qJD(5) * t23 + t309; t4, t35 * qJD(3), t39 * qJD(4) + t30 * qJD(5) + (-t341 + t387) * t294 + ((t297 + t308) * t367 + m(5) * t168 / 0.2e1 + (t297 + t313) * t369) * t371 + ((t343 - t388) * qJD(3) - t224 * t337) * t223 + t257, t39 * qJD(3), -t315 + t30 * qJD(3) + m(7) * t141 * qJD(5) + (t342 * qJD(5) * t224 + (t344 * qJD(5) - t337) * t222) * t225 (-t222 * t293 - t223 * t295) * m(7); qJD(4) * t11 + qJD(5) * t2 - qJD(6) * t19 - t263, -qJD(4) * t48 - qJD(5) * t29 - t257, qJD(4) * t55 + qJD(5) * t18 - qJD(6) * t64, t242, t157 * qJD(6) + t243 + (-mrSges(7,2) * t311 - Ifges(6,6) * t224 + t271 * t222 + t377 * t226 + t214) * qJD(5), qJD(5) * t157 - t256; -qJD(3) * t11 - qJD(5) * t8 + t223 * t290 + t310, t48 * qJD(3), -t242, 0, t377 * qJD(5) + t290 - t314 (qJD(1) * t223 + qJD(5)) * t222 * m(7); -qJD(3) * t2 + qJD(4) * t8 + qJD(6) * t24 - t262, qJD(3) * t29 + t315, -t243, t314, t212 * qJD(6), t254; qJD(3) * t19 - qJD(4) * t291 - qJD(5) * t24 - t309, 0, t256, -qJD(1) * t291, -t254, 0;];
Cq  = t13;

% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:59
% EndTime: 2019-12-31 21:04:10
% DurationCPUTime: 4.41s
% Computational Cost: add. (4215->468), mult. (9327->610), div. (0->0), fcn. (7119->4), ass. (0->213)
t250 = sin(qJ(2));
t249 = sin(qJ(3));
t236 = t249 * mrSges(6,1);
t251 = cos(qJ(3));
t323 = t251 * mrSges(6,2);
t366 = -t236 + t323;
t150 = t366 * t250;
t380 = -t150 / 0.2e1;
t370 = Ifges(5,4) + Ifges(4,5);
t235 = t250 * qJ(4);
t221 = t249 * t235;
t310 = t250 * t251;
t348 = pkin(3) + pkin(4);
t296 = t348 * t310;
t84 = -t221 - t296;
t307 = t348 * t249;
t309 = t251 * qJ(4);
t379 = t307 - t309;
t369 = Ifges(4,6) + Ifges(6,6);
t327 = Ifges(5,5) * t251;
t328 = Ifges(6,4) * t251;
t378 = t327 + t328 + (Ifges(6,2) + Ifges(5,3)) * t249;
t242 = Ifges(6,4) * t249;
t279 = Ifges(6,1) * t251 + t242;
t240 = Ifges(5,5) * t249;
t280 = Ifges(5,1) * t251 + t240;
t329 = Ifges(4,4) * t249;
t281 = Ifges(4,1) * t251 - t329;
t377 = t279 + t280 + t281;
t355 = m(6) / 0.2e1;
t376 = 0.2e1 * t355;
t375 = -t236 / 0.2e1;
t374 = -t250 / 0.2e1;
t373 = t250 / 0.2e1;
t252 = cos(qJ(2));
t372 = -t252 / 0.2e1;
t371 = t252 / 0.2e1;
t368 = pkin(7) - qJ(5);
t196 = mrSges(6,1) * t251 + mrSges(6,2) * t249;
t147 = t196 * t250;
t316 = qJ(4) * t249;
t275 = -pkin(3) * t251 - t316;
t192 = -pkin(2) + t275;
t284 = t251 * mrSges(5,1) + t249 * mrSges(5,3);
t367 = -m(5) * t192 + t284;
t365 = Ifges(5,6) * t249 + t251 * t370;
t244 = Ifges(4,4) * t251;
t278 = -Ifges(4,2) * t249 + t244;
t364 = Ifges(5,3) * t251 - t240;
t363 = Ifges(6,2) * t251 - t242;
t362 = -Ifges(4,1) * t249 - t244;
t361 = Ifges(5,2) + Ifges(4,3) + Ifges(6,3);
t319 = t252 * Ifges(6,5);
t107 = t250 * t279 + t319;
t109 = -t252 * Ifges(5,4) + t250 * t280;
t111 = -t252 * Ifges(4,5) + t250 * t281;
t298 = -Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t353 = -Ifges(6,5) / 0.2e1;
t290 = t353 - t298;
t360 = -t290 * t252 + t107 / 0.2e1 + t109 / 0.2e1 + t111 / 0.2e1;
t359 = 0.2e1 * m(6);
t358 = -m(5) / 0.2e1;
t357 = m(5) / 0.2e1;
t356 = -m(6) / 0.2e1;
t354 = -mrSges(5,3) / 0.2e1;
t352 = -Ifges(5,6) / 0.2e1;
t333 = pkin(7) * t252;
t210 = pkin(2) * t250 - t333;
t295 = -pkin(6) * t249 - pkin(3);
t35 = (-qJ(5) * t252 - t210) * t251 + (-pkin(4) + t295) * t250;
t351 = t35 / 0.2e1;
t312 = t249 * t250;
t220 = qJ(5) * t312;
t308 = t251 * t252;
t193 = -pkin(2) * t252 - pkin(7) * t250 - pkin(1);
t314 = t193 * t249;
t88 = pkin(6) * t308 + t314;
t54 = t220 + t88;
t350 = -t54 / 0.2e1;
t313 = t210 * t251;
t72 = t250 * t295 - t313;
t349 = t72 / 0.2e1;
t283 = mrSges(5,1) * t249 - mrSges(5,3) * t251;
t149 = t283 * t250;
t346 = t149 / 0.2e1;
t345 = -t284 / 0.2e1;
t344 = t196 / 0.2e1;
t198 = t368 * t251;
t343 = t198 / 0.2e1;
t342 = -t249 / 0.2e1;
t341 = t249 / 0.2e1;
t339 = -t251 / 0.2e1;
t335 = m(6) * t379;
t334 = m(6) * t249;
t332 = -mrSges(5,2) + mrSges(6,3);
t197 = pkin(3) * t249 - t309;
t265 = pkin(6) + t197;
t113 = t265 * t250;
t222 = mrSges(6,3) * t312;
t238 = t252 * mrSges(6,2);
t180 = t222 - t238;
t237 = t252 * mrSges(5,3);
t301 = mrSges(5,2) * t312;
t191 = -t237 - t301;
t68 = t314 + (pkin(6) * t251 - qJ(4)) * t252;
t42 = t220 + t68;
t264 = -pkin(6) - t379;
t65 = t264 * t250;
t9 = (-t180 - t191) * t252 + (-t149 + t150) * t310 + m(5) * (-t113 * t310 - t252 * t68) + m(6) * (-t252 * t42 + t310 * t65);
t325 = qJD(1) * t9;
t320 = t252 * mrSges(6,1);
t114 = t265 * t252;
t151 = t283 * t252;
t152 = t366 * t252;
t199 = mrSges(4,1) * t249 + mrSges(4,2) * t251;
t153 = t252 * t199;
t181 = mrSges(4,2) * t252 - mrSges(4,3) * t312;
t311 = t249 * t252;
t182 = mrSges(6,2) * t250 + mrSges(6,3) * t311;
t183 = -mrSges(4,2) * t250 - mrSges(4,3) * t311;
t184 = -mrSges(6,3) * t310 + t320;
t185 = -mrSges(4,1) * t252 - mrSges(4,3) * t310;
t186 = mrSges(5,1) * t252 + mrSges(5,2) * t310;
t223 = mrSges(6,3) * t308;
t187 = -t250 * mrSges(6,1) - t223;
t188 = mrSges(4,1) * t250 - mrSges(4,3) * t308;
t300 = mrSges(5,2) * t308;
t189 = -t250 * mrSges(5,1) + t300;
t190 = -mrSges(5,2) * t311 + mrSges(5,3) * t250;
t262 = Ifges(6,5) * t374 + t290 * t250 + t370 * t373 + t371 * t377;
t287 = Ifges(5,6) * t373 + t278 * t372 + t369 * t374 + t371 * t378;
t225 = Ifges(5,5) * t310;
t101 = -Ifges(5,6) * t252 + Ifges(5,3) * t312 + t225;
t226 = Ifges(6,4) * t310;
t103 = Ifges(6,2) * t312 + t252 * Ifges(6,6) + t226;
t105 = -t252 * Ifges(4,6) + t250 * t278;
t288 = t101 / 0.2e1 + t103 / 0.2e1 - t105 / 0.2e1;
t297 = Ifges(4,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t289 = t352 + t297;
t306 = pkin(6) * t311 - t193 * t251;
t246 = t252 * pkin(3);
t53 = qJ(5) * t310 - t306;
t34 = pkin(4) * t252 + t246 - t53;
t91 = -pkin(6) * t310 + t210 * t249;
t71 = t235 + t91;
t43 = qJ(5) * t311 + t71;
t66 = t264 * t252;
t69 = t246 + t306;
t90 = pkin(6) * t312 + t313;
t3 = t114 * t149 + t66 * t150 + t113 * t151 + t65 * t152 + t43 * t180 + t91 * t181 + t42 * t182 + t88 * t183 + t35 * t184 + t90 * t185 + t72 * t186 + t34 * t187 - t306 * t188 + t69 * t189 + t68 * t190 + t71 * t191 + m(4) * (-t306 * t90 + t88 * t91) + m(6) * (t34 * t35 + t42 * t43 + t65 * t66) + m(5) * (t113 * t114 + t68 * t71 + t69 * t72) + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t252 + t360 * t251 + (t252 * t289 + t288) * t249) * t252 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t250 + pkin(6) * t153 + t262 * t251 + (-t289 * t250 + t287) * t249 + (Ifges(3,1) - Ifges(3,2) + (m(4) * pkin(6) + t199) * pkin(6) - t361) * t252) * t250;
t318 = t3 * qJD(1);
t145 = pkin(3) * t310 + t221;
t146 = t284 * t250;
t285 = t251 * mrSges(4,1) - t249 * mrSges(4,2);
t148 = t285 * t250;
t154 = t364 * t250;
t155 = t363 * t250;
t202 = Ifges(4,2) * t251 + t329;
t156 = t250 * t202;
t157 = -Ifges(6,1) * t312 + t226;
t158 = -Ifges(5,1) * t312 + t225;
t159 = t362 * t250;
t224 = Ifges(5,6) * t310;
t4 = t113 * t146 - t65 * t147 + t145 * t149 + t84 * t150 + t53 * t180 + t54 * t184 + t224 * t372 + (-t185 + t186) * t88 - (t181 + t191) * t306 + m(6) * (t34 * t54 + t42 * t53 + t65 * t84) + m(5) * (t113 * t145 - t306 * t68 + t69 * t88) + (pkin(6) * t148 + (-t88 * mrSges(4,3) - t68 * mrSges(5,2) + t42 * mrSges(6,3) + t157 / 0.2e1 + t158 / 0.2e1 + t159 / 0.2e1 + t297 * t252 + t288) * t251 + (-t306 * mrSges(4,3) - t69 * mrSges(5,2) + t34 * mrSges(6,3) + t154 / 0.2e1 + t155 / 0.2e1 + t156 / 0.2e1 - t360) * t249) * t250;
t317 = t4 * qJD(1);
t15 = (t249 * t180 - t251 * t184 + m(6) * (t249 * t42 - t34 * t251)) * t250;
t315 = qJD(1) * t15;
t304 = qJD(2) * t249;
t303 = m(6) * t310;
t302 = t66 * t355;
t299 = -mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1;
t160 = t251 * t348 + pkin(2) + t316;
t194 = t368 * t249;
t254 = (t113 * t197 + t145 * t192) * t357 + (t160 * t84 - t379 * t65 + (t34 + t53) * t198 + (-t42 + t54) * t194) * t355 - pkin(2) * t148 / 0.2e1 + t145 * t345 - t160 * t147 / 0.2e1 + t379 * t380 + t192 * t146 / 0.2e1 - t194 * t180 / 0.2e1 + t197 * t346 + t184 * t343 + t65 * t375 + t84 * t344 - t365 * t252 / 0.4e1;
t203 = Ifges(6,1) * t249 - t328;
t204 = Ifges(5,1) * t249 - t327;
t247 = t249 ^ 2;
t248 = t251 ^ 2;
t255 = pkin(6) * t199 / 0.2e1 + (mrSges(6,3) * t343 + t240 / 0.4e1 + t242 / 0.4e1 - t202 / 0.4e1 - t363 / 0.4e1 - t364 / 0.4e1 + (Ifges(4,1) / 0.4e1 + Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1) * t251) * t251 + (t194 * mrSges(6,3) / 0.2e1 + t362 / 0.4e1 - t204 / 0.4e1 - t203 / 0.4e1 - t244 / 0.4e1 + (Ifges(4,2) / 0.4e1 + Ifges(6,2) / 0.4e1 + Ifges(5,3) / 0.4e1) * t249 + (-Ifges(4,4) / 0.4e1 + Ifges(6,4) / 0.4e1 + Ifges(5,5) / 0.4e1) * t251) * t249 + (mrSges(5,2) + mrSges(4,3)) * pkin(7) * (-t248 / 0.2e1 - t247 / 0.2e1);
t256 = (-pkin(3) * t72 + qJ(4) * t71) * t358 + (qJ(4) * t43 - t348 * t35) * t356 + pkin(3) * t189 / 0.2e1 + t348 * t187 / 0.2e1 + mrSges(6,1) * t351 - t43 * mrSges(6,2) / 0.2e1 + t71 * t354 + mrSges(5,1) * t349 - t90 * mrSges(4,1) / 0.2e1 + t91 * mrSges(4,2) / 0.2e1;
t257 = (-t53 / 0.2e1 - t34 / 0.2e1) * mrSges(6,3) + (-t306 / 0.2e1 + t69 / 0.2e1) * mrSges(5,2) + (t186 / 0.2e1 - t185 / 0.2e1 + (t69 - t306) * t357) * pkin(7) + t107 / 0.4e1 + t109 / 0.4e1 + t111 / 0.4e1 - t154 / 0.4e1 - t155 / 0.4e1 - t156 / 0.4e1 + t113 * t354 + t65 * mrSges(6,2) / 0.2e1;
t258 = (t350 + t42 / 0.2e1) * mrSges(6,3) + (t88 / 0.2e1 - t68 / 0.2e1) * mrSges(5,2) + (-t191 / 0.2e1 - t181 / 0.2e1 + (-t68 + t88) * t357) * pkin(7) + t101 / 0.4e1 + t103 / 0.4e1 - t105 / 0.4e1 + t157 / 0.4e1 + t158 / 0.4e1 + t159 / 0.4e1 + t113 * mrSges(5,1) / 0.2e1;
t2 = ((0.3e1 / 0.4e1 * Ifges(4,6) + t352 + 0.3e1 / 0.4e1 * Ifges(6,6)) * t252 + t258) * t249 + ((0.3e1 / 0.4e1 * Ifges(6,5) + t298) * t252 + t257) * t251 + (-Ifges(5,2) / 0.2e1 - Ifges(4,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t255) * t250 + t254 + t256 + (-t182 / 0.2e1 - t190 / 0.2e1) * qJ(4);
t5 = pkin(2) * t199 + t379 * t196 - t192 * t283 + t202 * t341 + t378 * t251 / 0.2e1 + t367 * t197 + (t335 - t366) * t160 + (t203 + t204 - t362 + t278) * t339 + (-t364 - t363 + t377) * t342;
t274 = t2 * qJD(1) - t5 * qJD(2);
t21 = (m(6) * t160 + t196 + t367) * t249;
t259 = (-t113 * t249 + (-t192 * t250 - t333) * t251) * t358 + (t160 * t310 - t198 * t252 + t249 * t65) * t356;
t269 = m(5) * t349 + m(6) * t351;
t6 = t300 - t223 + (t346 + t380) * t249 + ((t345 - t196 / 0.2e1) * t251 + t299) * t250 + t259 + t269;
t273 = -qJD(1) * t6 + qJD(2) * t21;
t261 = m(6) * ((-t194 * t250 - t42) * t251 + (t198 * t250 - t34) * t249);
t10 = (t238 / 0.2e1 + t180 / 0.2e1) * t251 + (-t320 / 0.2e1 + t184 / 0.2e1) * t249 + t302 - t261 / 0.2e1;
t29 = m(6) * (-t194 * t249 - t198 * t251) + (t247 + t248) * mrSges(6,3);
t272 = -qJD(1) * t10 + qJD(2) * t29;
t20 = t147 + (t296 / 0.4e1 + t221 / 0.4e1 - t84 / 0.4e1) * t359;
t26 = (-t379 / 0.4e1 + t309 / 0.4e1 - t307 / 0.4e1) * t359 + t366;
t271 = qJD(1) * t20 - qJD(2) * t26;
t263 = -0.2e1 * qJ(4) * t252 + t88;
t260 = -t237 - t238 + t263 * t357 + (t220 + t263) * t355;
t268 = m(6) * t350 + t358 * t88;
t14 = t260 + t268;
t211 = mrSges(6,2) + mrSges(5,3) + (m(5) + m(6)) * qJ(4);
t270 = qJD(1) * t14 + qJD(3) * t211;
t267 = -mrSges(5,2) * pkin(3) + mrSges(6,3) * t348 - Ifges(6,5);
t266 = qJ(4) * t332 - t369;
t161 = (qJD(1) * t310 + t304) * m(6);
t52 = m(6) * t198 + (m(5) * pkin(7) - t332) * t251;
t37 = t379 * t355 - t335 / 0.2e1;
t13 = t222 + t260 - t268 - t301;
t11 = t184 * t342 + t180 * t339 + t261 / 0.2e1 + t302 + (t323 / 0.2e1 + t375) * t252;
t7 = t149 * t342 + t150 * t341 + t299 * t250 - t259 + t269 + (t284 / 0.2e1 + t344) * t310;
t1 = ((Ifges(4,6) / 0.4e1 + Ifges(6,6) / 0.4e1) * t252 + t258) * t249 + (t319 / 0.4e1 + t257) * t251 + t255 * t250 + t254 - t256 + (t182 + t190) * qJ(4) / 0.2e1 + t361 * t373 + (Ifges(5,6) / 0.2e1 - t369 / 0.2e1) * t311 + (t353 + t370 / 0.2e1) * t308;
t8 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t9 + qJD(5) * t15, t318 + t1 * qJD(3) + t7 * qJD(4) + t11 * qJD(5) + (t72 * mrSges(5,2) - t90 * mrSges(4,3) - t35 * mrSges(6,3) + t262) * t304 + (-pkin(2) * t153 + t192 * t151 + t160 * t152 + t198 * t182 + t194 * t187 + t66 * t196 + (t160 * t66 + t194 * t35 + t198 * t43) * t376 + (t71 * mrSges(5,2) + t91 * mrSges(4,3) - t43 * mrSges(6,3) - t287) * t251 + ((t183 + t190) * t251 + (-t188 + t189) * t249 + m(5) * (t249 * t72 + t251 * t71) + m(4) * (-t249 * t90 + t251 * t91)) * pkin(7) + (Ifges(3,5) + (t203 / 0.2e1 + t204 / 0.2e1 - t362 / 0.2e1) * t251 + (-t364 / 0.2e1 - t363 / 0.2e1 - t202 / 0.2e1) * t249 + (-m(4) * pkin(2) - mrSges(3,1) - t285) * pkin(6)) * t252 + (pkin(6) * mrSges(3,2) + t251 * t289 - Ifges(3,6)) * t250 - t367 * t114) * qJD(2), t1 * qJD(2) + t13 * qJD(4) + t317 + (-t88 * mrSges(4,1) - t88 * mrSges(5,1) - t54 * mrSges(6,1) + t306 * mrSges(4,2) + t53 * mrSges(6,2) - t306 * mrSges(5,3) + t224 + 0.2e1 * (-pkin(3) * t88 - qJ(4) * t306) * t357 + (qJ(4) * t53 - t348 * t54) * t376 + (t266 * t251 + (-t267 - t370) * t249) * t250) * qJD(3), qJD(2) * t7 + qJD(3) * t13 + t325, qJD(2) * t11 + t315; qJD(3) * t2 - qJD(4) * t6 - qJD(5) * t10 - t318, -qJD(3) * t5 + qJD(4) * t21 + qJD(5) * t29, t52 * qJD(4) + t37 * qJD(5) + t274 + (m(6) * (-qJ(4) * t194 - t198 * t348) - t198 * mrSges(6,1) - t194 * mrSges(6,2) + t267 * t251 + t266 * t249 + (m(5) * t275 - t284 - t285) * pkin(7) + t365) * qJD(3), qJD(3) * t52 + t273, qJD(3) * t37 + t272; -qJD(2) * t2 + qJD(4) * t14 + qJD(5) * t20 - t317, -qJD(5) * t26 - t274, t211 * qJD(4), t270, t271; qJD(2) * t6 - qJD(3) * t14 - qJD(5) * t303 - t325, -qJD(5) * t334 - t273, -t270, 0, -t161; qJD(2) * t10 - qJD(3) * t20 + qJD(4) * t303 - t315, qJD(3) * t26 + qJD(4) * t334 - t272, -t271, t161, 0;];
Cq = t8;

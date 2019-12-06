% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP8
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:20
% EndTime: 2019-12-05 16:58:32
% DurationCPUTime: 4.76s
% Computational Cost: add. (4141->445), mult. (10786->628), div. (0->0), fcn. (9677->8), ass. (0->219)
t381 = m(5) / 0.2e1;
t222 = sin(qJ(3));
t225 = cos(qJ(3));
t220 = sin(pkin(5));
t223 = sin(qJ(2));
t306 = t220 * t223;
t319 = cos(pkin(5));
t159 = t222 * t306 - t225 * t319;
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t253 = pkin(4) * t224 + qJ(5) * t221;
t164 = t253 * t222;
t303 = t221 * t222;
t172 = mrSges(5,2) * t225 - mrSges(5,3) * t303;
t212 = t225 * mrSges(6,3);
t179 = -mrSges(6,2) * t303 - t212;
t275 = -t179 / 0.2e1 - t172 / 0.2e1;
t301 = t222 * t224;
t174 = -mrSges(5,1) * t225 - mrSges(5,3) * t301;
t175 = mrSges(6,1) * t225 + mrSges(6,2) * t301;
t276 = t174 / 0.2e1 - t175 / 0.2e1;
t300 = t224 * t225;
t346 = pkin(8) * t222;
t184 = -pkin(3) * t225 - pkin(2) - t346;
t304 = t221 * t184;
t122 = pkin(7) * t300 + t304;
t348 = pkin(7) * t224;
t271 = -qJ(5) + t348;
t96 = t225 * t271 + t304;
t321 = t122 - t96;
t302 = t221 * t225;
t308 = t184 * t224;
t121 = -pkin(7) * t302 + t308;
t286 = pkin(7) * t221 + pkin(4);
t97 = t225 * t286 - t308;
t322 = t121 + t97;
t380 = -m(6) / 0.2e1;
t160 = t222 * t319 + t225 * t306;
t226 = cos(qJ(2));
t299 = t224 * t226;
t76 = t160 * t221 + t220 * t299;
t305 = t220 * t226;
t288 = t221 * t305;
t77 = t160 * t224 - t288;
t415 = t276 * t77 + (t159 * t164 + t321 * t76 + t322 * t77) * t380 - t275 * t76;
t259 = t224 * mrSges(6,1) + t221 * mrSges(6,3);
t261 = t224 * mrSges(5,1) - t221 * mrSges(5,2);
t395 = t259 + t261;
t414 = -t395 / 0.2e1;
t216 = t221 ^ 2;
t218 = t224 ^ 2;
t413 = t216 + t218;
t182 = -pkin(3) - t253;
t120 = (t221 * t223 + t225 * t299) * t220;
t313 = t120 * t224;
t119 = -t224 * t306 + t225 * t288;
t315 = t119 * t221;
t238 = (t313 + t315) * pkin(8);
t289 = t222 * t305;
t403 = mrSges(6,2) + mrSges(5,3);
t412 = -t403 * (-t315 / 0.2e1 - t313 / 0.2e1) + (-pkin(3) * t289 + t238) * t381 - (t182 * t289 + t238) * t380;
t411 = t172 + t179;
t379 = m(6) / 0.2e1;
t407 = 0.2e1 * t379;
t406 = 0.2e1 * t381;
t405 = m(5) + m(6);
t404 = mrSges(5,1) + mrSges(6,1);
t402 = Ifges(6,4) + Ifges(5,5);
t401 = Ifges(5,6) - Ifges(6,6);
t399 = m(6) * t182 - t259;
t213 = Ifges(6,5) * t221;
t398 = Ifges(6,1) * t224 + t213;
t214 = Ifges(5,4) * t224;
t193 = t221 * Ifges(5,1) + t214;
t397 = -Ifges(5,2) * t221 + t214;
t340 = mrSges(6,3) * t224;
t342 = mrSges(6,1) * t221;
t258 = -t340 + t342;
t318 = qJ(5) * t224;
t350 = pkin(4) * t221;
t188 = -t318 + t350;
t351 = m(6) * t188;
t396 = t258 + t351;
t393 = -t413 * t346 / 0.2e1;
t388 = t224 * Ifges(6,3) - t213;
t392 = t398 - t388;
t344 = pkin(8) * t225;
t198 = pkin(3) * t222 - t344;
t307 = t198 * t224;
t125 = pkin(7) * t303 + t307;
t171 = t221 * t198;
t126 = -pkin(7) * t301 + t171;
t390 = -t125 * t221 + t126 * t224;
t104 = -t222 * t286 - t307;
t99 = -t222 * t271 + t171;
t389 = t104 * t221 + t224 * t99;
t326 = t225 * mrSges(4,2);
t189 = t222 * mrSges(4,1) + t326;
t337 = Ifges(6,5) * t224;
t192 = t221 * Ifges(6,1) - t337;
t386 = t397 + t192 + t193;
t211 = m(6) * qJ(5) + mrSges(6,3);
t138 = -t225 * Ifges(6,4) + t222 * t398;
t338 = Ifges(5,4) * t221;
t257 = Ifges(5,1) * t224 - t338;
t239 = t257 * t222;
t140 = -t225 * Ifges(5,5) + t239;
t292 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t385 = -t292 * t225 + t138 / 0.2e1 + t140 / 0.2e1;
t375 = mrSges(6,3) / 0.2e1;
t376 = -mrSges(5,2) / 0.2e1;
t384 = (-pkin(4) * t119 + qJ(5) * t120) * t379 + (t376 + t375) * t120;
t383 = 0.2e1 * m(6);
t382 = 0.2e1 * t159;
t378 = mrSges(5,1) / 0.2e1;
t377 = -mrSges(6,1) / 0.2e1;
t374 = -Ifges(6,5) / 0.2e1;
t373 = pkin(7) * mrSges(5,1);
t372 = pkin(7) * mrSges(5,2);
t370 = t122 / 0.2e1;
t368 = t159 / 0.2e1;
t178 = -mrSges(6,2) * t302 + mrSges(6,3) * t222;
t364 = t178 / 0.2e1;
t362 = -t182 / 0.2e1;
t361 = t259 / 0.2e1;
t360 = -t221 / 0.2e1;
t359 = t221 / 0.2e1;
t357 = t222 / 0.2e1;
t356 = -t224 / 0.2e1;
t355 = t224 / 0.2e1;
t341 = mrSges(5,2) * t224;
t343 = mrSges(5,1) * t221;
t260 = t341 + t343;
t166 = t260 * t222;
t349 = pkin(7) * t166;
t245 = pkin(7) + t188;
t142 = t245 * t222;
t335 = t142 * mrSges(6,1);
t334 = t142 * mrSges(6,3);
t330 = t222 * mrSges(6,1);
t329 = t222 * mrSges(4,2);
t325 = t225 * Ifges(5,6);
t324 = t225 * t77;
t323 = t77 * t224;
t320 = -t261 - mrSges(4,1);
t309 = t159 * t221;
t11 = t405 * (-t309 * t76 + (t160 - t323) * t159);
t316 = t11 * qJD(1);
t12 = m(4) * (t159 * t222 + t160 * t225 - t306) * t305 + t405 * (t119 * t76 + t77 * t120 + t159 * t289);
t314 = t12 * qJD(1);
t298 = t413 * pkin(8) * t159;
t296 = mrSges(6,2) * t300;
t295 = t104 * t379;
t291 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t290 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t280 = -t301 / 0.2e1;
t279 = t301 / 0.4e1;
t206 = Ifges(6,5) * t301;
t134 = -t225 * Ifges(6,6) + Ifges(6,3) * t303 + t206;
t136 = t222 * t397 - t325;
t278 = t134 / 0.2e1 - t136 / 0.2e1;
t272 = -mrSges(6,2) * qJ(5) - Ifges(5,6);
t266 = mrSges(6,2) * pkin(4) - t402;
t191 = t224 * Ifges(5,2) + t338;
t254 = Ifges(6,3) * t221 + t337;
t165 = t258 * t222;
t173 = -mrSges(5,2) * t222 - mrSges(5,3) * t302;
t176 = mrSges(5,1) * t222 - mrSges(5,3) * t300;
t177 = t296 - t330;
t227 = (t165 / 0.2e1 + t166 / 0.2e1) * t160 + (t177 / 0.2e1 - t176 / 0.2e1) * t76 + (t364 + t173 / 0.2e1) * t77 + (pkin(7) * t160 * t222 - t125 * t76 + t126 * t77) * t381 + (t104 * t76 + t142 * t160 + t77 * t99) * t379;
t167 = t258 * t225;
t168 = t260 * t225;
t231 = t167 / 0.2e1 + t168 / 0.2e1 + t275 * t224 + t276 * t221;
t233 = m(5) * (pkin(7) * t225 + t121 * t221 - t122 * t224);
t143 = t245 * t225;
t237 = m(6) * (-t221 * t97 - t224 * t96 + t143);
t1 = (t326 / 0.2e1 - t189 / 0.2e1 + (mrSges(4,1) / 0.2e1 + t261 / 0.2e1 + t361) * t222) * t305 + t231 * t159 + (t237 / 0.4e1 + t233 / 0.4e1) * t382 + t227 - t412;
t135 = Ifges(6,6) * t222 + t225 * t254;
t137 = Ifges(5,6) * t222 + t225 * t397;
t139 = Ifges(6,4) * t222 + t225 * t398;
t141 = Ifges(5,5) * t222 + t225 * t257;
t5 = -pkin(2) * t189 + t104 * t175 + t121 * t176 + t122 * t173 + t125 * t174 + t126 * t172 + t142 * t167 + t143 * t165 + t97 * t177 + t96 * t178 + t99 * t179 + m(6) * (t104 * t97 + t142 * t143 + t96 * t99) + m(5) * (t121 * t125 + t122 * t126) + (t349 + Ifges(4,4) * t225 + t385 * t224 + (-t290 * t225 + t278) * t221) * t225 + (-Ifges(4,4) * t222 + pkin(7) * t168 + (t139 / 0.2e1 + t141 / 0.2e1 + t292 * t222) * t224 + (t135 / 0.2e1 - t137 / 0.2e1 + t290 * t222) * t221 + (m(5) * pkin(7) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3)) * t225) * t222;
t252 = t1 * qJD(1) + t5 * qJD(2);
t6 = (-mrSges(5,1) / 0.2e1 + t377) * t119 + (t159 * t414 + t403 * (t76 * t359 + t323 / 0.2e1)) * t222 + t384 + t415;
t205 = Ifges(6,6) * t301;
t8 = t121 * t172 - t122 * t174 + t122 * t175 + t164 * t165 + m(6) * (t121 * t96 + t122 * t97 + t142 * t164) + t121 * t179 - t225 * t205 / 0.2e1 + ((t206 / 0.2e1 - t122 * mrSges(5,3) + t335 - t96 * mrSges(6,2) + t325 / 0.2e1 + (t373 - t214 / 0.2e1) * t222 + t278) * t224 + (t121 * mrSges(5,3) - t97 * mrSges(6,2) + t334 + (-t372 + (t374 + Ifges(5,4) / 0.2e1) * t221) * t222 + (-Ifges(6,1) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t301 - t385) * t221) * t222;
t251 = -t6 * qJD(1) + t8 * qJD(2);
t23 = t165 * t301 - m(6) * (-t142 * t301 - t225 * t96) + t225 * t179;
t24 = (t119 / 0.4e1 + t159 * t279 + t324 / 0.4e1) * t383;
t250 = qJD(1) * t24 + qJD(2) * t23;
t230 = (-t142 * t221 + (-t182 * t222 - t344) * t224) * t379 + t165 * t360;
t20 = t296 + (-t259 * t355 + t377) * t222 + t295 - t230;
t78 = t399 * t221;
t249 = qJD(2) * t20 + qJD(3) * t78;
t248 = t142 * t188 + t164 * t182;
t30 = -t212 + (t304 / 0.4e1 - t122 / 0.4e1 + (t348 / 0.4e1 - qJ(5) / 0.2e1) * t225) * t383;
t246 = qJD(2) * t30 + qJD(4) * t211;
t244 = Ifges(5,2) / 0.4e1 + Ifges(6,3) / 0.4e1 - Ifges(5,1) / 0.4e1 - Ifges(6,1) / 0.4e1;
t241 = t164 * t361 - t188 * t165 / 0.2e1;
t240 = t261 * t222;
t18 = -t188 * t259 + t191 * t360 + t254 * t356 - pkin(3) * t260 + t396 * t182 + (t257 + t392) * t359 + t386 * t355;
t229 = (-pkin(4) * t104 + qJ(5) * t99) * t379 - pkin(4) * t177 / 0.2e1 + qJ(5) * t364 + t104 * t377 + t125 * t378 + t126 * t376 + t99 * t375;
t4 = t248 * t380 + (t334 / 0.2e1 - t140 / 0.4e1 - t138 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(6,4) + 0.3e1 / 0.4e1 * Ifges(5,5)) * t225 + (-t121 / 0.2e1 - t97 / 0.2e1) * mrSges(6,2) + (t322 * t380 + t276) * pkin(8)) * t224 + (-t335 / 0.2e1 + t136 / 0.4e1 - t134 / 0.4e1 - t206 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(5,6)) * t225 + (-t122 / 0.2e1 + t96 / 0.2e1) * mrSges(6,2) + (t321 * t380 - t275) * pkin(8)) * t221 + ((pkin(3) * t378 + mrSges(6,1) * t362 - t372 / 0.2e1 - t213 / 0.4e1 + t191 / 0.4e1 + t388 / 0.4e1 + t244 * t224) * t224 + (pkin(3) * t376 + mrSges(6,3) * t362 + t193 / 0.4e1 + t192 / 0.4e1 + t214 / 0.4e1 - t373 / 0.2e1 - t244 * t221 + (0.3e1 / 0.4e1 * Ifges(5,4) + t374) * t224) * t221 + t291 + t403 * pkin(8) * (t218 / 0.2e1 + t216 / 0.2e1)) * t222 + t229 + t241;
t9 = (t350 / 0.4e1 - t318 / 0.4e1 - t188 / 0.4e1) * m(6) * t382;
t235 = t9 * qJD(1) + t4 * qJD(2) - t18 * qJD(3);
t234 = t259 * t357;
t219 = t225 ^ 2;
t217 = t222 ^ 2;
t185 = t217 * pkin(7) * t305;
t183 = (m(6) * pkin(8) + mrSges(6,2)) * t224;
t50 = m(6) * t309;
t28 = (t304 + (-0.2e1 * qJ(5) + t348) * t225) * t379 + m(6) * t370 + t179;
t25 = (-t159 * t301 + t119 - t324) * t379;
t21 = -t259 * t280 - t330 / 0.2e1 + t295 + t230;
t10 = (t341 / 0.2e1 + t343 / 0.2e1 + t342 / 0.2e1 + t351 / 0.2e1 - t340 / 0.2e1) * t159 + (t260 + t396) * t368;
t7 = t240 * t368 + t159 * t234 - t404 * t119 / 0.2e1 + t384 + t403 * (t77 * t280 - t76 * t303 / 0.2e1) - t415;
t3 = t191 * t280 + t182 * t234 - t241 + t349 / 0.2e1 + t229 + t248 * t379 - t388 * t301 / 0.4e1 + (t221 * t290 + t224 * t292) * t225 + t291 * t222 - t221 * t136 / 0.4e1 - pkin(3) * t240 / 0.2e1 + t142 * t258 / 0.2e1 + (-Ifges(6,1) * t303 + t134 + t206) * t221 / 0.4e1 - (-t401 * t221 + t402 * t224) * t225 / 0.4e1 + (t254 - t193) * t303 / 0.4e1 + t392 * t279 + t393 * mrSges(5,3) + (t239 + t140 + t138) * t224 / 0.4e1 - t386 * t303 / 0.4e1 + ((t221 * t321 + t224 * t322) * t379 - t276 * t224 + t411 * t360) * pkin(8) + ((t370 - t96 / 0.2e1) * t221 + t322 * t355 + t393) * mrSges(6,2);
t2 = (t237 / 0.2e1 + t233 / 0.2e1 + t231) * t159 + t227 + t289 * t414 - t189 * t305 + t412;
t13 = [qJD(2) * t12 + qJD(3) * t11, t2 * qJD(3) + t7 * qJD(4) + t25 * qJD(5) + t314 + (((-t225 * mrSges(4,1) - mrSges(3,1) + t329) * t223 + (-mrSges(3,2) + (t165 + t166) * t222 + (t217 + t219) * mrSges(4,3)) * t226) * t220 + t185 * t406 + t142 * t289 * t407 + m(4) * (t185 + (pkin(7) * t219 * t226 - pkin(2) * t223) * t220) + (t122 * t406 + t96 * t407 + t411) * t120 + (-t121 * t406 + t97 * t407 - t174 + t175) * t119) * qJD(2), t2 * qJD(2) + t10 * qJD(4) - t50 * qJD(5) + t316 + ((-t259 + t320) * t160 + (-pkin(3) * t160 - t298) * t406 + (t160 * t182 - t298) * t407 + (-t403 * t413 + mrSges(4,2)) * t159) * qJD(3), t7 * qJD(2) + t10 * qJD(3) + (-t404 * t77 + (mrSges(5,2) - mrSges(6,3)) * t76) * qJD(4) + ((-pkin(4) * t77 - qJ(5) * t76) * qJD(4) / 0.2e1 + t77 * qJD(5) / 0.2e1) * t383, m(6) * t77 * qJD(4) + t25 * qJD(2) - t50 * qJD(3); qJD(3) * t1 - qJD(4) * t6 - qJD(5) * t24 - t314, qJD(3) * t5 + qJD(4) * t8 - qJD(5) * t23, t3 * qJD(4) + t21 * qJD(5) + t252 + (t135 * t356 + t137 * t355 - Ifges(4,6) * t222 + t182 * t167 - pkin(3) * t168 + pkin(7) * t329 + ((t173 + t178) * t224 + (-t176 + t177) * t221 + m(6) * t389 + m(5) * t390) * pkin(8) + (Ifges(4,5) + (t192 / 0.2e1 + t193 / 0.2e1) * t224 + (-t388 / 0.2e1 - t191 / 0.2e1) * t221 + (-m(5) * pkin(3) + t320) * pkin(7)) * t225 + (t139 + t141) * t359 + (t402 * t221 + t401 * t224) * t357 + t399 * t143 + t390 * mrSges(5,3) + t389 * mrSges(6,2)) * qJD(3), t3 * qJD(3) + t28 * qJD(5) + t251 + (t205 + (t221 * t266 + t224 * t272) * t222 + (-m(6) * pkin(4) - t404) * t122 + (-mrSges(5,2) + t211) * t121) * qJD(4), qJD(3) * t21 + qJD(4) * t28 - t250; -qJD(2) * t1 - qJD(4) * t9 - t316, -qJD(4) * t4 - qJD(5) * t20 - t252, qJD(4) * t18 - qJD(5) * t78, t183 * qJD(5) - t235 + (-t266 * t224 + (Ifges(6,6) + t272) * t221 + (-m(6) * t253 - t395) * pkin(8)) * qJD(4), qJD(4) * t183 - t249; qJD(2) * t6 + qJD(3) * t9, qJD(3) * t4 + qJD(5) * t30 - t251, t235, t211 * qJD(5), t246; qJD(2) * t24, qJD(3) * t20 - qJD(4) * t30 + t250, t249, -t246, 0;];
Cq = t13;

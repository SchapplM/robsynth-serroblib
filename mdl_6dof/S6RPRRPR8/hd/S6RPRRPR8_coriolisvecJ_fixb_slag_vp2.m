% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 16:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:20:10
% EndTime: 2018-11-23 16:20:22
% DurationCPUTime: 11.96s
% Computational Cost: add. (9789->635), mult. (22107->895), div. (0->0), fcn. (14794->8), ass. (0->294)
t253 = sin(qJ(3));
t256 = cos(qJ(3));
t282 = pkin(3) * t256 + pkin(8) * t253;
t228 = t282 * qJD(1);
t257 = -pkin(1) - pkin(7);
t239 = qJD(1) * t257 + qJD(2);
t255 = cos(qJ(4));
t252 = sin(qJ(4));
t316 = t252 * t256;
t160 = t228 * t255 - t239 * t316;
t341 = -qJ(5) - pkin(8);
t284 = qJD(4) * t341;
t430 = -(qJ(5) * t253 * t255 + pkin(4) * t256) * qJD(1) - t160 - qJD(5) * t252 + t255 * t284;
t314 = t255 * t256;
t161 = t228 * t252 + t239 * t314;
t248 = t253 * qJD(1);
t294 = t252 * t248;
t302 = qJD(5) * t255;
t429 = qJ(5) * t294 - t252 * t284 + t161 - t302;
t249 = sin(pkin(10));
t250 = cos(pkin(10));
t219 = t249 * t255 + t250 * t252;
t202 = t219 * qJD(1);
t177 = t253 * t202;
t201 = t219 * qJD(4);
t313 = t177 + t201;
t268 = t249 * t252 - t250 * t255;
t394 = t268 * t253;
t178 = qJD(1) * t394;
t389 = qJD(4) * t268;
t312 = t178 + t389;
t412 = t249 * t429 + t250 * t430;
t411 = t249 * t430 - t250 * t429;
t251 = sin(qJ(6));
t254 = cos(qJ(6));
t244 = t248 + qJD(4);
t308 = qJD(3) * t255;
t310 = qJD(1) * t256;
t224 = -t252 * t310 + t308;
t225 = qJD(3) * t252 + t255 * t310;
t155 = t224 * t249 + t225 * t250;
t417 = pkin(9) * t155;
t230 = pkin(3) * t253 - pkin(8) * t256 + qJ(2);
t207 = t230 * qJD(1);
t229 = t253 * t239;
t210 = qJD(3) * pkin(8) + t229;
t146 = t207 * t252 + t210 * t255;
t111 = qJ(5) * t224 + t146;
t105 = t249 * t111;
t145 = t207 * t255 - t210 * t252;
t110 = -qJ(5) * t225 + t145;
t95 = pkin(4) * t244 + t110;
t52 = t250 * t95 - t105;
t34 = pkin(5) * t244 - t417 + t52;
t283 = t224 * t250 - t225 * t249;
t404 = pkin(9) * t283;
t317 = t250 * t111;
t53 = t249 * t95 + t317;
t38 = t53 + t404;
t10 = t251 * t34 + t254 * t38;
t291 = t253 * t308;
t301 = qJD(3) * qJD(4);
t303 = qJD(4) * t256;
t167 = t255 * t301 + (-t252 * t303 - t291) * qJD(1);
t289 = t255 * t303;
t309 = qJD(3) * t253;
t408 = -t252 * t309 + t289;
t168 = -qJD(1) * t408 - t252 * t301;
t103 = -t167 * t249 + t168 * t250;
t104 = t167 * t250 + t168 * t249;
t407 = -t155 * t251 + t254 * t283;
t28 = qJD(6) * t407 + t103 * t251 + t104 * t254;
t307 = qJD(3) * t256;
t286 = qJD(1) * t307;
t85 = t155 * t254 + t251 * t283;
t29 = -qJD(6) * t85 + t103 * t254 - t104 * t251;
t300 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t286;
t351 = Ifges(7,4) * t85;
t237 = qJD(6) + t244;
t359 = -t237 / 0.2e1;
t79 = Ifges(7,4) * t407;
t37 = Ifges(7,1) * t85 + Ifges(7,5) * t237 + t79;
t379 = -t85 / 0.2e1;
t381 = -t407 / 0.2e1;
t217 = qJD(3) * t282 + qJD(2);
t192 = t217 * qJD(1);
t262 = -qJD(4) * t146 + t192 * t255;
t48 = -qJ(5) * t167 - qJD(5) * t225 + (pkin(4) * qJD(1) - t239 * t252) * t307 + t262;
t293 = t239 * t307;
t304 = qJD(4) * t255;
t305 = qJD(4) * t252;
t76 = t192 * t252 + t207 * t304 - t210 * t305 + t255 * t293;
t55 = qJ(5) * t168 + qJD(5) * t224 + t76;
t15 = -t249 * t55 + t250 * t48;
t11 = pkin(5) * t286 - pkin(9) * t104 + t15;
t16 = t249 * t48 + t250 * t55;
t12 = pkin(9) * t103 + t16;
t9 = -t251 * t38 + t254 * t34;
t2 = qJD(6) * t9 + t11 * t251 + t12 * t254;
t3 = -qJD(6) * t10 + t11 * t254 - t12 * t251;
t403 = mrSges(7,1) * t3 - t2 * mrSges(7,2);
t318 = t239 * t256;
t211 = -qJD(3) * pkin(3) - t318;
t158 = -pkin(4) * t224 + qJD(5) + t211;
t92 = -pkin(5) * t283 + t158;
t428 = t300 + t403 + (Ifges(7,5) * t407 - Ifges(7,6) * t85) * t359 + (t10 * t85 + t407 * t9) * mrSges(7,3) + (-Ifges(7,2) * t85 + t37 + t79) * t381 - t92 * (mrSges(7,1) * t85 + mrSges(7,2) * t407) + (Ifges(7,1) * t407 - t351) * t379;
t427 = -pkin(5) * t310 + pkin(9) * t312 + t412;
t426 = pkin(9) * t313 - t411;
t36 = Ifges(7,2) * t407 + Ifges(7,6) * t237 + t351;
t424 = t36 / 0.2e1;
t422 = Ifges(5,3) + Ifges(6,3);
t188 = t219 * t256;
t393 = qJD(1) * t268 - qJD(3) * t188 + t253 * t389;
t190 = t268 * t256;
t392 = -qJD(3) * t190 - t201 * t253 - t202;
t337 = Ifges(5,4) * t225;
t141 = Ifges(5,2) * t224 + Ifges(5,6) * t244 + t337;
t213 = Ifges(5,4) * t224;
t142 = Ifges(5,1) * t225 + Ifges(5,5) * t244 + t213;
t272 = t145 * t255 + t146 * t252;
t335 = Ifges(5,4) * t255;
t277 = -Ifges(5,2) * t252 + t335;
t336 = Ifges(5,4) * t252;
t279 = Ifges(5,1) * t255 - t336;
t280 = mrSges(5,1) * t252 + mrSges(5,2) * t255;
t333 = Ifges(5,6) * t252;
t334 = Ifges(5,5) * t255;
t353 = t255 / 0.2e1;
t355 = -t252 / 0.2e1;
t356 = t244 / 0.2e1;
t360 = t225 / 0.2e1;
t418 = -t272 * mrSges(5,3) + t141 * t355 + t142 * t353 + t211 * t280 + (-t333 + t334) * t356 + t277 * t224 / 0.2e1 + t279 * t360;
t416 = qJD(1) / 0.2e1;
t234 = t341 * t252;
t235 = t341 * t255;
t163 = t234 * t250 + t235 * t249;
t134 = -pkin(9) * t219 + t163;
t164 = t234 * t249 - t235 * t250;
t135 = -pkin(9) * t268 + t164;
t64 = t134 * t251 + t135 * t254;
t415 = -qJD(6) * t64 + t251 * t426 + t254 * t427;
t63 = t134 * t254 - t135 * t251;
t414 = qJD(6) * t63 + t251 * t427 - t254 * t426;
t413 = t155 * Ifges(6,4);
t183 = -pkin(4) * t294 + t229;
t410 = pkin(4) * t305 + pkin(5) * t313 - t183;
t321 = Ifges(4,5) * qJD(3);
t339 = Ifges(4,4) * t253;
t409 = t321 / 0.2e1 + (Ifges(4,1) * t256 - t339) * t416 + t418;
t74 = Ifges(6,2) * t283 + Ifges(6,6) * t244 + t413;
t406 = t74 / 0.2e1;
t371 = -t283 / 0.2e1;
t287 = -Ifges(4,6) * qJD(3) / 0.2e1;
t56 = -t103 * mrSges(6,1) + mrSges(6,2) * t104;
t8 = -t29 * mrSges(7,1) + mrSges(7,2) * t28;
t405 = -t56 - t8;
t402 = Ifges(6,4) * t283;
t246 = pkin(4) * t250 + pkin(5);
t349 = pkin(4) * t249;
t194 = t246 * t251 + t254 * t349;
t58 = -t110 * t249 - t317;
t40 = t58 - t404;
t59 = t110 * t250 - t105;
t41 = t59 - t417;
t401 = -qJD(6) * t194 + t251 * t41 - t254 * t40;
t193 = t246 * t254 - t251 * t349;
t400 = qJD(6) * t193 - t251 * t40 - t254 * t41;
t187 = t219 * t253;
t121 = -t187 * t254 + t251 * t394;
t399 = qJD(6) * t121 + t251 * t393 + t254 * t392;
t123 = -t187 * t251 - t254 * t394;
t398 = -qJD(6) * t123 - t251 * t392 + t254 * t393;
t397 = qJ(2) * (m(3) + m(4));
t391 = Ifges(5,5) * t167 + Ifges(5,6) * t168;
t315 = t253 * t257;
t174 = t230 * t252 + t255 * t315;
t390 = Ifges(6,5) * t104 + Ifges(6,6) * t103 + t286 * t422 + t391;
t77 = -t252 * t293 + t262;
t273 = -t252 * t77 + t255 * t76;
t171 = -mrSges(5,2) * t244 + mrSges(5,3) * t224;
t172 = mrSges(5,1) * t244 - mrSges(5,3) * t225;
t269 = -t171 * t252 - t172 * t255;
t388 = -m(5) * t272 + t269;
t299 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t338 = Ifges(4,4) * t256;
t387 = -t299 * t244 - t145 * mrSges(5,1) - t52 * mrSges(6,1) - t9 * mrSges(7,1) - t225 * Ifges(5,5) - t224 * Ifges(5,6) - t287 + (-Ifges(4,2) * t253 + t338) * t416 - t237 * Ifges(7,3) - t85 * Ifges(7,5) - t407 * Ifges(7,6) - t155 * Ifges(6,5) - t283 * Ifges(6,6) + t10 * mrSges(7,2) + t146 * mrSges(5,2) + t53 * mrSges(6,2) - t422 * t356;
t386 = t77 * mrSges(5,1) + t15 * mrSges(6,1) - t76 * mrSges(5,2) - t16 * mrSges(6,2);
t384 = m(5) / 0.2e1;
t383 = t28 / 0.2e1;
t382 = t29 / 0.2e1;
t380 = t407 / 0.2e1;
t378 = t85 / 0.2e1;
t375 = t103 / 0.2e1;
t374 = t104 / 0.2e1;
t122 = -t188 * t254 + t190 * t251;
t373 = t122 / 0.2e1;
t124 = -t188 * t251 - t190 * t254;
t372 = t124 / 0.2e1;
t370 = t283 / 0.2e1;
t369 = -t155 / 0.2e1;
t368 = t155 / 0.2e1;
t367 = t167 / 0.2e1;
t366 = t168 / 0.2e1;
t365 = -t188 / 0.2e1;
t364 = -t190 / 0.2e1;
t363 = -t224 / 0.2e1;
t361 = -t225 / 0.2e1;
t358 = t237 / 0.2e1;
t357 = -t244 / 0.2e1;
t352 = m(6) * t158;
t350 = pkin(4) * t225;
t261 = -qJD(4) * t174 + t217 * t255;
t285 = -t252 * t257 + pkin(4);
t71 = qJ(5) * t291 + (qJ(5) * t305 + qJD(3) * t285 - t302) * t256 + t261;
t306 = qJD(3) * t257;
t290 = t256 * t306;
t295 = t217 * t252 + t230 * t304 + t255 * t290;
t78 = -qJ(5) * t289 + (-qJD(5) * t256 + (qJ(5) * qJD(3) - qJD(4) * t257) * t253) * t252 + t295;
t33 = t249 * t71 + t250 * t78;
t340 = mrSges(4,1) * t253;
t332 = qJ(2) * mrSges(4,1);
t331 = qJ(2) * mrSges(4,2);
t115 = t177 * t254 - t178 * t251;
t149 = t219 * t254 - t251 * t268;
t81 = -qJD(6) * t149 - t201 * t254 + t251 * t389;
t323 = t115 - t81;
t116 = t177 * t251 + t178 * t254;
t148 = -t219 * t251 - t254 * t268;
t80 = qJD(6) * t148 - t201 * t251 - t254 * t389;
t322 = t116 - t80;
t319 = qJD(3) * mrSges(4,2);
t216 = t255 * t230;
t147 = -qJ(5) * t314 + t253 * t285 + t216;
t159 = -qJ(5) * t316 + t174;
t87 = t147 * t249 + t159 * t250;
t311 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t224 - mrSges(5,2) * t225 - mrSges(4,3) * t310;
t298 = t252 * t315;
t247 = -pkin(4) * t255 - pkin(3);
t226 = t239 * t309;
t288 = -t321 / 0.2e1;
t138 = -pkin(4) * t168 + t226;
t32 = -t249 * t78 + t250 * t71;
t86 = t147 * t250 - t159 * t249;
t220 = pkin(4) * t316 - t256 * t257;
t281 = mrSges(5,1) * t255 - mrSges(5,2) * t252;
t278 = Ifges(5,1) * t252 + t335;
t276 = Ifges(5,2) * t255 + t336;
t274 = Ifges(5,5) * t252 + Ifges(5,6) * t255;
t61 = pkin(5) * t253 + pkin(9) * t190 + t86;
t62 = -pkin(9) * t188 + t87;
t30 = -t251 * t62 + t254 * t61;
t31 = t251 * t61 + t254 * t62;
t271 = t145 * t252 - t146 * t255;
t150 = mrSges(5,1) * t286 - mrSges(5,3) * t167;
t151 = -mrSges(5,2) * t286 + mrSges(5,3) * t168;
t270 = -t150 * t252 + t151 * t255;
t169 = pkin(4) * t408 + t253 * t306;
t232 = -mrSges(4,3) * t248 - t319;
t227 = (mrSges(4,2) * t256 + t340) * qJD(1);
t179 = pkin(5) * t268 + t247;
t173 = t216 - t298;
t157 = pkin(5) * t188 + t220;
t131 = qJD(3) * t394 - t201 * t256;
t129 = t219 * t309 + t256 * t389;
t120 = mrSges(6,1) * t244 - mrSges(6,3) * t155;
t119 = -mrSges(6,2) * t244 + mrSges(6,3) * t283;
t114 = pkin(5) * t155 + t350;
t113 = -t252 * t290 + t261;
t112 = -qJD(4) * t298 + t295;
t109 = -mrSges(5,1) * t168 + mrSges(5,2) * t167;
t97 = Ifges(5,1) * t167 + Ifges(5,4) * t168 + Ifges(5,5) * t286;
t96 = Ifges(5,4) * t167 + Ifges(5,2) * t168 + Ifges(5,6) * t286;
t91 = mrSges(6,1) * t286 - mrSges(6,3) * t104;
t90 = -mrSges(6,2) * t286 + mrSges(6,3) * t103;
t89 = -pkin(5) * t129 + t169;
t88 = -mrSges(6,1) * t283 + mrSges(6,2) * t155;
t75 = Ifges(6,1) * t155 + Ifges(6,5) * t244 + t402;
t68 = mrSges(7,1) * t237 - mrSges(7,3) * t85;
t67 = -mrSges(7,2) * t237 + mrSges(7,3) * t407;
t60 = -pkin(5) * t103 + t138;
t50 = Ifges(6,1) * t104 + Ifges(6,4) * t103 + Ifges(6,5) * t286;
t49 = Ifges(6,4) * t104 + Ifges(6,2) * t103 + Ifges(6,6) * t286;
t45 = -qJD(6) * t124 + t129 * t254 - t131 * t251;
t43 = qJD(6) * t122 + t129 * t251 + t131 * t254;
t39 = -mrSges(7,1) * t407 + mrSges(7,2) * t85;
t24 = -mrSges(7,2) * t286 + mrSges(7,3) * t29;
t23 = mrSges(7,1) * t286 - mrSges(7,3) * t28;
t22 = pkin(9) * t129 + t33;
t21 = pkin(5) * t307 - pkin(9) * t131 + t32;
t7 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t286;
t6 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t286;
t5 = -qJD(6) * t31 + t21 * t254 - t22 * t251;
t4 = qJD(6) * t30 + t21 * t251 + t22 * t254;
t1 = [(Ifges(7,1) * t124 + Ifges(7,4) * t122) * t383 + (t227 + ((2 * mrSges(3,3)) + t340 + 0.2e1 * t397) * qJD(1)) * qJD(2) + (Ifges(7,4) * t124 + Ifges(7,2) * t122) * t382 + (Ifges(6,6) * t375 + Ifges(7,6) * t382 + Ifges(7,5) * t383 + Ifges(6,5) * t374 + (0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1) + Ifges(7,3) / 0.2e1 + t299) * t286 + t386 + t403) * t253 + (t129 * t53 - t131 * t52 + t15 * t190 - t16 * t188) * mrSges(6,3) + (-Ifges(6,1) * t190 - Ifges(6,4) * t188) * t374 + (-Ifges(6,4) * t190 - Ifges(6,2) * t188) * t375 + t138 * (mrSges(6,1) * t188 - mrSges(6,2) * t190) + (t279 * t367 + t277 * t366 + t96 * t355 + t97 * t353 - t257 * t109 + qJD(1) * qJD(2) * mrSges(4,2) + (-t252 * t76 - t255 * t77) * mrSges(5,3) + (t211 * t281 + t276 * t363 + t278 * t361 + t274 * t357 + t142 * t355 - t255 * t141 / 0.2e1 + t271 * mrSges(5,3)) * qJD(4) + ((0.2e1 * t332 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t334 / 0.2e1 - t333 / 0.2e1) * t256 + Ifges(6,5) * t364 + Ifges(6,6) * t365 + Ifges(7,5) * t372 + Ifges(7,6) * t373) * qJD(1) + (-m(5) * t257 + t280) * t229 + t287 + t257 * t232 - t387) * qJD(3)) * t256 + (t10 * t45 + t122 * t2 - t124 * t3 - t43 * t9) * mrSges(7,3) + (t300 + t390 + t391) * t253 / 0.2e1 + t129 * t406 + m(5) * (t112 * t146 + t113 * t145 + t173 * t77 + t174 * t76) + (Ifges(7,1) * t43 + Ifges(7,4) * t45) * t378 + (Ifges(7,4) * t43 + Ifges(7,2) * t45) * t380 + (Ifges(6,1) * t131 + Ifges(6,4) * t129) * t368 + (Ifges(6,4) * t131 + Ifges(6,2) * t129) * t370 + t7 * t372 + t6 * t373 + (Ifges(7,5) * t43 + Ifges(7,6) * t45) * t358 + t50 * t364 + t49 * t365 + t45 * t424 + t30 * t23 + t31 * t24 + (Ifges(6,5) * t131 + Ifges(6,6) * t129) * t356 + t43 * t37 / 0.2e1 + t4 * t67 + t5 * t68 + t89 * t39 + t87 * t90 + t86 * t91 + t92 * (-mrSges(7,1) * t45 + mrSges(7,2) * t43) + t33 * t119 + t32 * t120 + t60 * (-mrSges(7,1) * t122 + mrSges(7,2) * t124) + t131 * t75 / 0.2e1 + t157 * t8 + t158 * (-mrSges(6,1) * t129 + mrSges(6,2) * t131) + t169 * t88 + t112 * t171 + t113 * t172 + m(7) * (t10 * t4 + t157 * t60 + t2 * t31 + t3 * t30 + t5 * t9 + t89 * t92) + m(6) * (t138 * t220 + t15 * t86 + t158 * t169 + t16 * t87 + t32 * t52 + t33 * t53) + t173 * t150 + (t288 + (-0.2e1 * t331 + 0.3e1 / 0.2e1 * t339) * qJD(1) + (m(5) * t211 - t311) * t257 - t409) * t309 + t174 * t151 + t220 * t56; t121 * t23 + t123 * t24 - t187 * t91 - t394 * t90 + t398 * t68 + t399 * t67 + t393 * t120 + t392 * t119 + (-t109 + (t171 * t255 - t172 * t252 + t232) * qJD(3) + t405) * t256 + (t269 * qJD(4) + (t39 + t88 - t311) * qJD(3) + t270) * t253 - m(5) * t271 * t307 + 0.2e1 * ((-qJD(4) * t272 + t273) * t384 + (t352 / 0.2e1 + (t211 - t318) * t384) * qJD(3)) * t253 + (t10 * t399 + t121 * t3 + t123 * t2 - t256 * t60 + t309 * t92 + t398 * t9) * m(7) + (-t138 * t256 - t15 * t187 - t16 * t394 + t392 * t53 + t393 * t52) * m(6) + (-t227 + (-mrSges(3,3) - t397) * qJD(1) + t388) * qJD(1); (t80 / 0.2e1 - t116 / 0.2e1) * t37 + (-t177 / 0.2e1 - t201 / 0.2e1) * t74 + (t81 / 0.2e1 - t115 / 0.2e1) * t36 + t273 * mrSges(5,3) + t270 * pkin(8) + (-Ifges(6,4) * t389 - Ifges(6,2) * t201) * t370 + (-Ifges(6,1) * t389 - Ifges(6,4) * t201) * t368 + (-t178 / 0.2e1 - t389 / 0.2e1) * t75 + (-Ifges(6,5) * t389 - Ifges(6,6) * t201) * t356 + (Ifges(6,1) * t219 - Ifges(6,4) * t268) * t374 + (Ifges(6,4) * t219 - Ifges(6,2) * t268) * t375 + (-t15 * t219 - t16 * t268 + t312 * t52 - t313 * t53) * mrSges(6,3) + t138 * (mrSges(6,1) * t268 + mrSges(6,2) * t219) - t268 * t49 / 0.2e1 + (Ifges(7,1) * t80 + Ifges(7,4) * t81) * t378 + (Ifges(7,1) * t116 + Ifges(7,4) * t115) * t379 + (Ifges(7,4) * t80 + Ifges(7,2) * t81) * t380 + (Ifges(7,4) * t116 + Ifges(7,2) * t115) * t381 + (Ifges(7,4) * t149 + Ifges(7,2) * t148) * t382 + (Ifges(7,1) * t149 + Ifges(7,4) * t148) * t383 + (Ifges(6,1) * t178 + Ifges(6,4) * t177) * t369 + (Ifges(6,4) * t178 + Ifges(6,2) * t177) * t371 + (Ifges(7,5) * t80 + Ifges(7,6) * t81) * t358 + (Ifges(7,5) * t116 + Ifges(7,6) * t115) * t359 + t276 * t366 + t278 * t367 + ((t88 + t352) * t252 * pkin(4) + t388 * pkin(8) + t418) * qJD(4) + t96 * t353 + (Ifges(6,5) * t178 + Ifges(6,6) * t177) * t357 + t63 * t23 + t64 * t24 - pkin(3) * t109 + t148 * t6 / 0.2e1 + t149 * t7 / 0.2e1 + t60 * (-mrSges(7,1) * t148 + mrSges(7,2) * t149) + (mrSges(6,1) * t313 - mrSges(6,2) * t312) * t158 + t163 * t91 + t164 * t90 - t161 * t171 - t160 * t172 + ((t288 + (t331 - t339 / 0.2e1) * qJD(1) + t409) * t253 + ((t338 / 0.2e1 - t332 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t253) * qJD(1) + t287 + t387) * t256 + (Ifges(6,5) * t219 + Ifges(7,5) * t149 - Ifges(6,6) * t268 + Ifges(7,6) * t148 + t274) * t307 / 0.2e1) * qJD(1) + t410 * t39 + t179 * t8 - t183 * t88 + ((-t232 - t319) * t256 + ((-mrSges(4,1) - t281) * qJD(3) + t311) * t253) * t239 + (mrSges(7,1) * t323 - mrSges(7,2) * t322) * t92 + (-t10 * t323 + t148 * t2 - t149 * t3 + t322 * t9) * mrSges(7,3) + t411 * t119 + t412 * t120 + (t138 * t247 + t15 * t163 - t158 * t183 + t16 * t164 + t411 * t53 + t412 * t52) * m(6) + t219 * t50 / 0.2e1 + t414 * t67 + t415 * t68 + (t10 * t414 + t179 * t60 + t2 * t64 + t3 * t63 + t410 * t92 + t415 * t9) * m(7) + t247 * t56 + t252 * t97 / 0.2e1 + (-pkin(3) * t226 + pkin(8) * t273 - t145 * t160 - t146 * t161 - t211 * t229) * m(5); t386 + t400 * t67 + (t10 * t400 - t114 * t92 + t193 * t3 + t194 * t2 + t401 * t9) * m(7) + t401 * t68 + ((t15 * t250 + t16 * t249) * pkin(4) - t158 * t350 - t52 * t58 - t53 * t59) * m(6) + t85 * t424 + (t145 * t224 + t146 * t225) * mrSges(5,3) + (-t225 * t88 + t249 * t90 + t250 * t91) * pkin(4) + t428 + (Ifges(5,5) * t224 + Ifges(6,5) * t283 - Ifges(5,6) * t225 - Ifges(6,6) * t155) * t357 + (t155 * t53 + t283 * t52) * mrSges(6,3) - t158 * (mrSges(6,1) * t155 + mrSges(6,2) * t283) + t155 * t406 + (-Ifges(6,2) * t155 + t402 + t75) * t371 + t390 + t141 * t360 + (Ifges(5,1) * t224 - t337) * t361 + (-Ifges(5,2) * t225 + t142 + t213) * t363 - t114 * t39 - t59 * t119 - t58 * t120 - t145 * t171 + t146 * t172 + t193 * t23 + t194 * t24 + (Ifges(6,1) * t283 - t413) * t369 - t211 * (mrSges(5,1) * t225 + mrSges(5,2) * t224); -t283 * t119 + t155 * t120 - t407 * t67 + t85 * t68 + (-t10 * t407 + t85 * t9 + t60) * m(7) + (t155 * t52 - t283 * t53 + t138) * m(6) - t405; t10 * t68 + t36 * t378 - t9 * t67 + t428;];
tauc  = t1(:);

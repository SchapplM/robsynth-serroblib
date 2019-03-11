% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:14
% EndTime: 2019-03-08 23:48:48
% DurationCPUTime: 16.24s
% Computational Cost: add. (8583->693), mult. (23288->930), div. (0->0), fcn. (18094->12), ass. (0->326)
t233 = cos(qJ(3));
t223 = sin(pkin(7));
t318 = qJD(2) * t223;
t297 = t233 * t318;
t434 = qJD(4) - t297;
t228 = sin(qJ(4));
t232 = cos(qJ(4));
t225 = cos(pkin(7));
t317 = qJD(2) * t225;
t291 = qJD(3) + t317;
t229 = sin(qJ(3));
t298 = t229 * t318;
t166 = t228 * t291 + t232 * t298;
t363 = t166 / 0.2e1;
t160 = qJD(6) + t166;
t367 = t160 / 0.2e1;
t260 = t232 * t291;
t165 = t228 * t298 - t260;
t227 = sin(qJ(6));
t231 = cos(qJ(6));
t125 = t165 * t227 + t231 * t434;
t371 = t125 / 0.2e1;
t124 = t165 * t231 - t227 * t434;
t373 = t124 / 0.2e1;
t427 = Ifges(7,5) * t371 + Ifges(7,6) * t373 + Ifges(7,3) * t367;
t366 = -t165 / 0.2e1;
t429 = Ifges(5,4) * t366;
t433 = Ifges(5,1) * t363 + t427 + t429;
t226 = cos(pkin(6));
t319 = qJD(1) * t226;
t301 = t223 * t319;
t230 = sin(qJ(2));
t224 = sin(pkin(6));
t320 = qJD(1) * t224;
t300 = t230 * t320;
t193 = pkin(9) * t318 + t300;
t234 = cos(qJ(2));
t201 = qJD(2) * pkin(2) + t234 * t320;
t335 = t225 * t229;
t396 = t233 * t193 + t201 * t335;
t115 = t229 * t301 + t396;
t289 = t228 * pkin(4) * t297 + t115;
t316 = qJD(4) * t228;
t292 = pkin(4) * t316 - qJD(5) * t228;
t339 = qJ(5) * t232;
t431 = -t289 + t292 + t434 * (pkin(11) * t228 - t339);
t376 = pkin(5) + pkin(10);
t214 = t376 * t232;
t377 = pkin(4) + pkin(11);
t184 = t229 * t193;
t114 = t233 * (t201 * t225 + t301) - t184;
t278 = pkin(3) * t229 - pkin(10) * t233;
t179 = t278 * t318;
t82 = -t228 * t114 + t179 * t232;
t430 = qJD(4) * t214 - (pkin(5) * t232 * t233 - t229 * t377) * t318 + t82;
t361 = t434 / 0.2e1;
t326 = t233 * t234;
t330 = t229 * t230;
t255 = -t225 * t330 + t326;
t151 = t255 * t320;
t288 = t223 * t300;
t121 = t151 * t228 - t232 * t288;
t337 = t223 * t233;
t191 = pkin(2) * t335 + pkin(9) * t337;
t174 = pkin(10) * t225 + t191;
t279 = -pkin(3) * t233 - pkin(10) * t229;
t175 = (-pkin(2) + t279) * t223;
t254 = t278 * qJD(3);
t180 = t223 * t254;
t338 = t223 * t229;
t217 = pkin(9) * t338;
t334 = t225 * t233;
t190 = pkin(2) * t334 - t217;
t181 = t190 * qJD(3);
t315 = qJD(4) * t232;
t59 = -t174 * t315 - t175 * t316 + t180 * t232 - t228 * t181;
t428 = -t121 - t59;
t426 = Ifges(6,5) / 0.2e1;
t293 = -qJ(5) * t228 - pkin(3);
t192 = -t232 * t377 + t293;
t213 = t376 * t228;
t141 = t192 * t231 + t213 * t227;
t425 = -qJD(6) * t141 - t431 * t227 + t430 * t231;
t140 = -t192 * t227 + t213 * t231;
t424 = qJD(6) * t140 + t430 * t227 + t431 * t231;
t304 = t228 * t338;
t284 = qJD(4) * t304;
t295 = qJD(3) * t337;
t144 = -t225 * t315 - t232 * t295 + t284;
t314 = t229 * qJD(3);
t296 = t223 * t314;
t261 = t377 * t296;
t423 = -pkin(5) * t144 - t261 + t428;
t189 = t225 * t228 + t232 * t338;
t145 = qJD(4) * t189 + t228 * t295;
t328 = t230 * t233;
t329 = t229 * t234;
t257 = t225 * t328 + t329;
t150 = t257 * t320;
t182 = t191 * qJD(3);
t243 = qJ(5) * t144 - qJD(5) * t189 + t182;
t422 = -t145 * t377 + t150 - t243;
t331 = t228 * t233;
t83 = t232 * t114 + t228 * t179;
t421 = -(-pkin(5) * t331 + qJ(5) * t229) * t318 - t83 - t376 * t316;
t98 = -pkin(3) * t291 - t114;
t241 = -t166 * qJ(5) + t98;
t47 = t165 * pkin(4) + t241;
t299 = t225 * t319;
t244 = (qJD(2) * t279 - t201) * t223 + t299;
t99 = pkin(10) * t291 + t115;
t45 = t228 * t99 - t232 * t244;
t259 = pkin(5) * t166 + t45;
t417 = qJD(5) + t259;
t27 = -t377 * t434 + t417;
t38 = t165 * t377 + t241;
t7 = -t227 * t38 + t231 * t27;
t8 = t227 * t27 + t231 * t38;
t420 = t7 * mrSges(7,1) + t98 * mrSges(5,2) - t8 * mrSges(7,2) - t47 * mrSges(6,3) + Ifges(5,5) * t361 + t433;
t354 = -qJD(3) / 0.2e1;
t249 = t255 * qJD(2);
t285 = t226 * t295;
t71 = (t201 * t334 - t184) * qJD(3) + (t224 * t249 + t285) * qJD(1);
t419 = qJD(4) * t244 + t71;
t46 = t228 * t244 + t232 * t99;
t40 = -qJ(5) * t434 - t46;
t418 = t40 * mrSges(6,1) - t46 * mrSges(5,3);
t399 = -qJD(5) - t45;
t39 = -pkin(4) * t434 - t399;
t158 = Ifges(6,6) * t165;
t92 = Ifges(6,4) * t434 - t166 * Ifges(6,2) + t158;
t416 = t92 / 0.2e1 - t45 * mrSges(5,3) - t39 * mrSges(6,1);
t415 = -t98 * mrSges(5,1) + t47 * mrSges(6,2) + Ifges(5,6) * t361 - t434 * t426 + (Ifges(5,2) + Ifges(6,3)) * t366 + (Ifges(5,4) + Ifges(6,6)) * t363;
t364 = -t166 / 0.2e1;
t365 = t165 / 0.2e1;
t402 = Ifges(6,4) - Ifges(5,5);
t414 = -Ifges(6,2) * t364 - Ifges(6,6) * t365 - t402 * t361 + t420 + t433;
t294 = qJD(3) * t318;
t282 = t233 * t294;
t131 = qJD(2) * t284 - qJD(4) * t260 - t232 * t282;
t139 = (t254 + t300) * t318;
t17 = t139 * t232 - t228 * t419 - t99 * t315;
t10 = -pkin(5) * t131 - qJD(2) * t261 - t17;
t132 = qJD(4) * t166 + t228 * t282;
t250 = t257 * qJD(2);
t286 = t226 * t296;
t72 = t396 * qJD(3) + (t224 * t250 + t286) * qJD(1);
t237 = qJ(5) * t131 - qJD(5) * t166 + t72;
t18 = t132 * t377 + t237;
t1 = qJD(6) * t7 + t10 * t227 + t18 * t231;
t413 = t1 * mrSges(7,2);
t2 = -qJD(6) * t8 + t10 * t231 - t18 * t227;
t412 = t2 * mrSges(7,1);
t370 = -t131 / 0.2e1;
t411 = -t132 / 0.2e1;
t410 = t298 / 0.2e1;
t409 = Ifges(4,6) * t354;
t116 = -t228 * t174 + t175 * t232;
t104 = pkin(4) * t337 - t116;
t70 = pkin(5) * t189 + pkin(11) * t337 + t104;
t188 = -t232 * t225 + t304;
t173 = t217 + (-pkin(2) * t233 - pkin(3)) * t225;
t248 = -qJ(5) * t189 + t173;
t81 = t188 * t377 + t248;
t25 = -t227 * t81 + t231 * t70;
t408 = qJD(6) * t25 + t227 * t423 - t231 * t422;
t26 = t227 * t70 + t231 * t81;
t407 = -qJD(6) * t26 + t227 * t422 + t231 * t423;
t401 = Ifges(6,5) - Ifges(5,6);
t398 = t150 - t182;
t397 = -t151 + t181;
t395 = t225 * t326 - t330;
t394 = t1 * t227 + t2 * t231;
t309 = Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t393 = -t309 * t434 - t416 + t420 + t427;
t251 = t174 * t316 - t175 * t315 - t228 * t180 - t232 * t181;
t41 = -t223 * (qJ(5) * t314 - qJD(5) * t233) + t251;
t157 = -t201 * t223 + t299;
t308 = t426 - Ifges(5,6) / 0.2e1;
t310 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t348 = Ifges(4,6) * t225;
t391 = t308 * t165 - t309 * t166 + t310 * t434 + t157 * mrSges(4,1) + t39 * mrSges(6,2) + t409 - (t348 + (Ifges(4,4) * t229 + Ifges(4,2) * t233) * t223) * qJD(2) / 0.2e1 + Ifges(5,5) * t363 + Ifges(5,6) * t366 + Ifges(6,4) * t364 + Ifges(6,5) * t365 - t115 * mrSges(4,3) - t40 * mrSges(6,3) - t45 * mrSges(5,1) - t46 * mrSges(5,2) + (Ifges(5,3) + Ifges(6,1)) * t361;
t351 = Ifges(7,4) * t227;
t269 = Ifges(7,2) * t231 + t351;
t350 = Ifges(7,4) * t231;
t271 = Ifges(7,1) * t227 + t350;
t274 = mrSges(7,1) * t231 - mrSges(7,2) * t227;
t276 = t227 * t7 - t231 * t8;
t356 = t165 * pkin(5);
t30 = -t40 - t356;
t347 = Ifges(7,6) * t231;
t349 = Ifges(7,5) * t227;
t359 = -t231 / 0.2e1;
t360 = -t227 / 0.2e1;
t368 = -t160 / 0.2e1;
t372 = -t125 / 0.2e1;
t374 = -t124 / 0.2e1;
t352 = Ifges(7,4) * t125;
t43 = Ifges(7,2) * t124 + Ifges(7,6) * t160 + t352;
t123 = Ifges(7,4) * t124;
t44 = Ifges(7,1) * t125 + Ifges(7,5) * t160 + t123;
t389 = mrSges(7,3) * t276 + (t347 + t349) * t368 + t269 * t374 + t271 * t372 + t30 * t274 + t359 * t43 + t360 * t44;
t388 = -Ifges(5,4) * t363 - Ifges(5,2) * t366 + Ifges(6,6) * t364 + Ifges(6,3) * t365 + t361 * t401 - t415;
t283 = t229 * t294;
t55 = qJD(6) * t124 + t132 * t227 + t231 * t283;
t56 = -qJD(6) * t125 + t132 * t231 - t227 * t283;
t15 = t55 * Ifges(7,1) + t56 * Ifges(7,4) - t131 * Ifges(7,5);
t386 = t15 / 0.2e1;
t385 = t43 / 0.2e1;
t384 = t44 / 0.2e1;
t383 = t55 / 0.2e1;
t382 = t56 / 0.2e1;
t381 = -Ifges(6,4) * t283 / 0.2e1 + Ifges(6,2) * t370 + Ifges(6,6) * t411;
t12 = -pkin(4) * t283 - t17;
t375 = m(6) * t12;
t52 = Ifges(7,5) * t55;
t51 = Ifges(7,6) * t56;
t353 = qJD(3) / 0.2e1;
t142 = -t224 * t395 - t226 * t337;
t344 = t142 * t72;
t136 = mrSges(6,1) * t165 - mrSges(6,3) * t434;
t64 = -mrSges(7,1) * t124 + mrSges(7,2) * t125;
t342 = -t136 + t64;
t340 = qJ(5) * t165;
t336 = t224 * t230;
t333 = t227 * t228;
t332 = t228 * t231;
t327 = t231 * t233;
t108 = -t131 * mrSges(6,1) + mrSges(6,2) * t283;
t109 = mrSges(5,1) * t283 + mrSges(5,3) * t131;
t325 = t108 - t109;
t107 = mrSges(6,1) * t132 - mrSges(6,3) * t283;
t110 = -mrSges(5,2) * t283 - mrSges(5,3) * t132;
t324 = t110 - t107;
t323 = -mrSges(4,1) * t291 + mrSges(5,1) * t165 + mrSges(5,2) * t166 + mrSges(4,3) * t298;
t134 = -mrSges(5,2) * t434 - mrSges(5,3) * t165;
t322 = t134 - t136;
t135 = mrSges(5,1) * t434 - mrSges(5,3) * t166;
t137 = mrSges(6,1) * t166 + mrSges(6,2) * t434;
t321 = t135 - t137;
t117 = t232 * t174 + t228 * t175;
t13 = -Ifges(7,3) * t131 + t51 + t52;
t311 = -Ifges(5,1) / 0.2e1 - Ifges(6,2) / 0.2e1;
t307 = -Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t306 = Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1;
t305 = t64 + t322;
t113 = -mrSges(6,2) * t165 - mrSges(6,3) * t166;
t302 = t113 + t323;
t215 = Ifges(4,4) * t297;
t290 = t157 * mrSges(4,2) + Ifges(4,1) * t410 + Ifges(4,5) * t291 / 0.2e1 + t215 / 0.2e1;
t287 = t318 * t336;
t281 = t412 - t413;
t280 = -t72 * mrSges(4,1) - t71 * mrSges(4,2);
t275 = -mrSges(4,1) * t233 + mrSges(4,2) * t229;
t273 = mrSges(7,1) * t227 + mrSges(7,2) * t231;
t272 = Ifges(7,1) * t231 - t351;
t270 = -Ifges(7,2) * t227 + t350;
t268 = Ifges(7,5) * t231 - Ifges(7,6) * t227;
t33 = -mrSges(7,1) * t131 - mrSges(7,3) * t55;
t34 = mrSges(7,2) * t131 + mrSges(7,3) * t56;
t265 = t227 * t34 + t231 * t33;
t87 = -mrSges(7,2) * t160 + mrSges(7,3) * t124;
t88 = mrSges(7,1) * t160 - mrSges(7,3) * t125;
t264 = -t227 * t88 + t231 * t87;
t263 = t228 * t40 + t232 * t39;
t262 = -t228 * t46 + t232 * t45;
t256 = t225 * t329 + t328;
t143 = t224 * t256 + t226 * t338;
t187 = -t223 * t224 * t234 + t225 * t226;
t101 = t143 * t228 - t187 * t232;
t53 = t101 * t231 - t142 * t227;
t54 = t101 * t227 + t142 * t231;
t102 = t143 * t232 + t187 * t228;
t103 = qJ(5) * t337 - t117;
t258 = -t188 * t227 + t223 * t327;
t146 = t188 * t231 + t227 * t337;
t16 = t228 * t139 + t232 * t419 - t316 * t99;
t11 = -qJ(5) * t283 - qJD(5) * t434 - t16;
t247 = -t17 * mrSges(5,1) + t16 * mrSges(5,2) - t12 * mrSges(6,2) + t11 * mrSges(6,3);
t239 = t306 * t166 - t308 * t434 + t415 - t418;
t210 = Ifges(6,1) * t283;
t209 = Ifges(4,5) * t282;
t208 = Ifges(5,3) * t283;
t205 = -pkin(4) * t232 + t293;
t186 = -qJ(5) * t315 + t292;
t178 = t275 * t318;
t177 = -mrSges(4,2) * t291 + mrSges(4,3) * t297;
t171 = (mrSges(4,1) * t229 + mrSges(4,2) * t233) * t294;
t156 = (t227 * t331 + t229 * t231) * t318;
t155 = (-t227 * t229 + t228 * t327) * t318;
t130 = Ifges(6,4) * t131;
t129 = Ifges(5,5) * t131;
t128 = Ifges(6,5) * t132;
t127 = Ifges(5,6) * t132;
t122 = t151 * t232 + t228 * t288;
t111 = pkin(4) * t166 + t340;
t100 = pkin(4) * t188 + t248;
t97 = t285 + (qJD(3) * t395 + t249) * t224;
t96 = t286 + (qJD(3) * t256 + t250) * t224;
t86 = -t297 * t339 + t289;
t85 = t166 * t377 + t340;
t84 = -pkin(5) * t188 - t103;
t78 = qJD(6) * t258 + t145 * t231 - t227 * t296;
t77 = qJD(6) * t146 + t145 * t227 + t231 * t296;
t74 = mrSges(5,1) * t132 - mrSges(5,2) * t131;
t73 = -mrSges(6,2) * t132 + mrSges(6,3) * t131;
t69 = -pkin(4) * t298 - t82;
t65 = -qJ(5) * t298 - t83;
t63 = -t131 * Ifges(5,1) - t132 * Ifges(5,4) + Ifges(5,5) * t283;
t62 = -t131 * Ifges(5,4) - t132 * Ifges(5,2) + Ifges(5,6) * t283;
t60 = Ifges(6,5) * t283 + t131 * Ifges(6,6) + t132 * Ifges(6,3);
t50 = pkin(4) * t145 + t243;
t48 = -pkin(4) * t296 - t59;
t37 = t46 - t356;
t32 = t228 * t287 - t143 * t316 + (qJD(4) * t187 + t97) * t232;
t31 = qJD(4) * t102 + t228 * t97 - t232 * t287;
t29 = -pkin(5) * t145 - t41;
t24 = pkin(4) * t132 + t237;
t21 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t20 = t227 * t37 + t231 * t85;
t19 = -t227 * t85 + t231 * t37;
t14 = Ifges(7,4) * t55 + Ifges(7,2) * t56 - Ifges(7,6) * t131;
t9 = -pkin(5) * t132 - t11;
t6 = qJD(6) * t53 + t227 * t31 + t231 * t96;
t5 = -qJD(6) * t54 - t227 * t96 + t231 * t31;
t3 = [t187 * t171 + t97 * t177 + t53 * t33 + t54 * t34 + t5 * t88 + t6 * t87 - t321 * t31 + (-mrSges(3,1) * t230 - mrSges(3,2) * t234) * qJD(2) ^ 2 * t224 + (t73 + t74) * t142 + t325 * t101 + t302 * t96 + t305 * t32 + (t21 + t324) * t102 + (t178 * t336 + (t142 * t233 - t143 * t229) * qJD(3) * mrSges(4,3)) * t318 + m(5) * (-t101 * t17 + t102 * t16 + t31 * t45 + t32 * t46 + t96 * t98 + t344) + m(6) * (t101 * t12 - t102 * t11 + t142 * t24 + t31 * t39 - t32 * t40 + t47 * t96) + m(7) * (t1 * t54 + t102 * t9 + t2 * t53 + t30 * t32 + t5 * t7 + t6 * t8) + m(4) * (-t114 * t96 + t115 * t97 + t344 + t143 * t71 + (qJD(1) * t187 + t157) * t287); (t388 + t418) * t145 + (-t414 + t416) * t144 + (t1 * t146 + t2 * t258 - t7 * t77 + t78 * t8) * mrSges(7,3) + (-Ifges(7,4) * t258 + Ifges(7,2) * t146 + Ifges(7,6) * t189) * t382 + (-Ifges(7,1) * t258 + Ifges(7,4) * t146 + Ifges(7,5) * t189) * t383 + t9 * (-mrSges(7,1) * t146 - mrSges(7,2) * t258) - t258 * t386 - t251 * t134 + (-Ifges(5,4) * t188 - Ifges(7,5) * t258 + Ifges(7,6) * t146 + (Ifges(5,1) + Ifges(7,3)) * t189) * t370 + (-t16 * t188 - t17 * t189) * mrSges(5,3) + (t11 * t188 + t12 * t189) * mrSges(6,1) + t132 * (-Ifges(6,6) * t189 + Ifges(6,3) * t188) / 0.2e1 + (t209 / 0.2e1 + t280) * t225 - t305 * t122 - t302 * t150 + (t100 * t24 + t103 * t11 + t104 * t12 + (-t150 + t50) * t47 + (t122 + t41) * t40 + (-t121 + t48) * t39) * m(6) + (Ifges(7,4) * t77 + Ifges(7,2) * t78) * t373 + t131 * (-Ifges(6,2) * t189 + Ifges(6,6) * t188) / 0.2e1 + (t63 + t13) * t189 / 0.2e1 + (Ifges(7,5) * t77 + Ifges(7,6) * t78) * t367 + (Ifges(7,1) * t77 + Ifges(7,4) * t78) * t371 + t321 * t121 + t323 * t182 + (t116 * t17 + t117 * t16 + t173 * t72 - t398 * t98 + (-t122 - t251) * t46 + t428 * t45) * m(5) + t189 * t381 + t77 * t384 + t78 * t385 + t25 * t33 + t26 * t34 + t397 * t177 + (t114 * t398 + t115 * t397 - t190 * t72 + t191 * t71) * m(4) + t29 * t64 + t30 * (-mrSges(7,1) * t78 + mrSges(7,2) * t77) + t84 * t21 + t100 * t73 + t103 * t107 + t104 * t108 + t50 * t113 + t116 * t109 + t117 * t110 + t407 * t88 + t408 * t87 + (t1 * t26 + t2 * t25 + t84 * t9 + t408 * t8 + t407 * t7 + (-t122 + t29) * t30) * m(7) + (t72 * mrSges(4,3) * t229 - pkin(2) * t171 + (t129 / 0.2e1 + t127 / 0.2e1 - t208 / 0.2e1 - t210 / 0.2e1 - t130 / 0.2e1 - t128 / 0.2e1 + t71 * mrSges(4,3) - t308 * t132 - t309 * t131 + t247) * t233 + (-m(4) * t157 - t178 + (-m(4) * pkin(2) + t275) * t318) * t300 + (((t317 + t353) * Ifges(4,5) + (-qJD(2) * t190 - t114) * mrSges(4,3) + t290) * t233 + (t409 + (-0.3e1 / 0.2e1 * t348 - t191 * mrSges(4,3) - t309 * t189 + t308 * t188 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - t310) * t337) * qJD(2) + t391) * t229) * qJD(3) + (0.3e1 / 0.2e1 * t233 ^ 2 - 0.3e1 / 0.2e1 * t229 ^ 2) * Ifges(4,4) * t294) * t223 - t189 * t413 + t59 * t135 + t41 * t136 + t48 * t137 + t146 * t14 / 0.2e1 + (Ifges(5,4) * t189 - Ifges(5,2) * t188) * t411 + t189 * t412 + t173 * t74 + t188 * t60 / 0.2e1 - t188 * t62 / 0.2e1 + t72 * (mrSges(5,1) * t188 + mrSges(5,2) * t189) + t24 * (-mrSges(6,2) * t188 - mrSges(6,3) * t189); (t30 * (-mrSges(7,1) * t332 + mrSges(7,2) * t333) + (Ifges(7,4) * t333 + Ifges(7,2) * t332) * t373 + (Ifges(7,1) * t333 + Ifges(7,4) * t332) * t371 + (Ifges(7,5) * t333 + Ifges(7,6) * t332) * t367 + t333 * t384 + t332 * t385 + t262 * mrSges(5,3) + t263 * mrSges(6,1) + (m(5) * t262 + m(6) * t263) * pkin(10) + (t332 * t8 - t333 * t7) * mrSges(7,3) + (-t322 * pkin(10) + t388) * t228 + (-t321 * pkin(10) - t92 / 0.2e1 + t414) * t232) * qJD(4) + (Ifges(7,5) * t156 + Ifges(7,6) * t155) * t368 + (-t8 * t155 + t7 * t156) * mrSges(7,3) + (t9 * t274 - t55 * t271 / 0.2e1 - t56 * t269 / 0.2e1 + t15 * t360 + t14 * t359 - t72 * mrSges(5,1) + t24 * mrSges(6,2) - t11 * mrSges(6,1) + t16 * mrSges(5,3) + t62 / 0.2e1 - t60 / 0.2e1 + t307 * t132 + (t349 / 0.2e1 + t347 / 0.2e1 - t306) * t131 + (-t1 * t231 + t2 * t227) * mrSges(7,3) + (m(5) * t16 - m(6) * t11 + t324) * pkin(10) + (-t30 * t273 + t227 * t385 + t44 * t359 + t270 * t374 + t272 * t372 + t268 * t368 + (t227 * t8 + t231 * t7) * mrSges(7,3)) * qJD(6)) * t232 + t209 + t421 * t64 - m(5) * (t115 * t98 - t45 * t82 + t46 * t83) - m(6) * (t39 * t69 + t40 * t65 + t47 * t86) + (-m(5) * t72 - t74) * pkin(3) + m(6) * (t186 * t47 + t205 * t24) + (t12 * mrSges(6,1) - t17 * mrSges(5,3) + t13 / 0.2e1 + t381 + t63 / 0.2e1 - t24 * mrSges(6,3) + t72 * mrSges(5,2) + t52 / 0.2e1 + t51 / 0.2e1 - t306 * t132 + (-Ifges(7,3) / 0.2e1 + t311) * t131 + (-m(5) * t17 + t325 + t375) * pkin(10) + t281) * t228 + (Ifges(7,1) * t156 + Ifges(7,4) * t155) * t372 + (Ifges(7,4) * t156 + Ifges(7,2) * t155) * t374 + (-t86 + t186) * t113 + t424 * t87 + t425 * t88 + (t1 * t141 + t140 * t2 + t214 * t9 + t30 * t421 + t424 * t8 + t425 * t7) * m(7) - t323 * t115 + t280 + (((t354 + t317 / 0.2e1) * Ifges(4,6) + Ifges(4,4) * t410 + (-t228 * t402 - t232 * t401) * t353 - t391) * t229 + (-t215 / 0.2e1 + (-t225 * Ifges(4,5) / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t338) * qJD(2) + t114 * mrSges(4,3) + Ifges(4,5) * t354 + (t165 * t307 + t239) * t228 + (t165 * t306 + t166 * t311 - t393) * t232 - t290) * t233) * t318 - t83 * t134 - t82 * t135 - t65 * t136 - t69 * t137 + t140 * t33 + t141 * t34 - t155 * t43 / 0.2e1 - t156 * t44 / 0.2e1 - t30 * (-mrSges(7,1) * t155 + mrSges(7,2) * t156) - t114 * t177 + t205 * t73 + t214 * t21; -t394 * mrSges(7,3) + (qJ(5) * t9 - t19 * t7 - t20 * t8 + t30 * t417) * m(7) + t259 * t64 + t130 + t128 - t129 - t127 + t9 * t273 + (t21 - t107) * qJ(5) + t208 + t210 - t247 + (t239 + (t307 - t311) * t165 + t389) * t166 + t389 * qJD(6) + t342 * qJD(5) + t268 * t370 + t321 * t46 + t322 * t45 + t270 * t382 + t272 * t383 + t231 * t386 - (m(7) * t394 + t265 + (-m(7) * t276 + t264) * qJD(6)) * t377 + (-t158 / 0.2e1 + t429 + t393) * t165 + (-pkin(4) * t12 - qJ(5) * t11 - t111 * t47 - t39 * t46 + t399 * t40) * m(6) + t14 * t360 - t20 * t87 - t19 * t88 - pkin(4) * t108 - t111 * t113; -t342 * t434 + t264 * qJD(6) + (t113 + t264) * t166 + t375 - m(6) * (-t166 * t47 - t40 * t434) + t265 + t108 + (-t160 * t276 - t30 * t434 + t394) * m(7); -t30 * (mrSges(7,1) * t125 + mrSges(7,2) * t124) + (Ifges(7,1) * t124 - t352) * t372 + t43 * t371 + (Ifges(7,5) * t124 - Ifges(7,6) * t125) * t368 - t7 * t87 + t8 * t88 + (t124 * t7 + t125 * t8) * mrSges(7,3) + t281 + t13 + (-Ifges(7,2) * t125 + t123 + t44) * t374;];
tauc  = t3(:);

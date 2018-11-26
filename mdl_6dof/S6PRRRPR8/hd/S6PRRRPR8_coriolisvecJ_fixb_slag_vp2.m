% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:26:49
% EndTime: 2018-11-23 15:27:04
% DurationCPUTime: 15.92s
% Computational Cost: add. (8583->691), mult. (23288->927), div. (0->0), fcn. (18094->12), ass. (0->328)
t442 = -Ifges(6,5) / 0.2e1;
t234 = cos(qJ(3));
t224 = sin(pkin(7));
t321 = qJD(2) * t224;
t302 = t234 * t321;
t441 = qJD(4) - t302;
t229 = sin(qJ(4));
t233 = cos(qJ(4));
t226 = cos(pkin(7));
t320 = qJD(2) * t226;
t296 = qJD(3) + t320;
t230 = sin(qJ(3));
t303 = t230 * t321;
t165 = t229 * t296 + t233 * t303;
t364 = t165 / 0.2e1;
t159 = qJD(6) + t165;
t368 = t159 / 0.2e1;
t264 = t233 * t296;
t164 = t229 * t303 - t264;
t228 = sin(qJ(6));
t232 = cos(qJ(6));
t125 = t164 * t228 + t232 * t441;
t372 = t125 / 0.2e1;
t124 = t164 * t232 - t228 * t441;
t374 = t124 / 0.2e1;
t434 = Ifges(7,5) * t372 + Ifges(7,6) * t374 + Ifges(7,3) * t368;
t367 = -t164 / 0.2e1;
t436 = Ifges(5,4) * t367;
t440 = Ifges(5,1) * t364 + t434 + t436;
t227 = cos(pkin(6));
t322 = qJD(1) * t227;
t306 = t224 * t322;
t231 = sin(qJ(2));
t225 = sin(pkin(6));
t323 = qJD(1) * t225;
t305 = t231 * t323;
t193 = pkin(9) * t321 + t305;
t235 = cos(qJ(2));
t202 = qJD(2) * pkin(2) + t235 * t323;
t336 = t226 * t230;
t399 = t234 * t193 + t202 * t336;
t115 = t230 * t306 + t399;
t294 = t229 * pkin(4) * t302 + t115;
t318 = qJD(4) * t229;
t297 = pkin(4) * t318 - qJD(5) * t229;
t340 = qJ(5) * t233;
t438 = -t294 + t297 + t441 * (pkin(11) * t229 - t340);
t377 = pkin(5) + pkin(10);
t215 = t377 * t233;
t378 = pkin(4) + pkin(11);
t184 = t230 * t193;
t114 = t234 * (t202 * t226 + t306) - t184;
t282 = pkin(3) * t230 - pkin(10) * t234;
t179 = t282 * t321;
t82 = -t229 * t114 + t179 * t233;
t437 = qJD(4) * t215 - (pkin(5) * t233 * t234 - t230 * t378) * t321 + t82;
t362 = t441 / 0.2e1;
t327 = t234 * t235;
t331 = t230 * t231;
t259 = -t226 * t331 + t327;
t150 = t259 * t323;
t293 = t224 * t305;
t121 = t150 * t229 - t233 * t293;
t338 = t224 * t234;
t191 = pkin(2) * t336 + pkin(9) * t338;
t174 = pkin(10) * t226 + t191;
t283 = -pkin(3) * t234 - pkin(10) * t230;
t175 = (-pkin(2) + t283) * t224;
t258 = t282 * qJD(3);
t180 = t224 * t258;
t339 = t224 * t230;
t218 = pkin(9) * t339;
t335 = t226 * t234;
t190 = pkin(2) * t335 - t218;
t181 = t190 * qJD(3);
t317 = qJD(4) * t233;
t59 = -t174 * t317 - t175 * t318 + t180 * t233 - t229 * t181;
t435 = -t121 - t59;
t298 = -qJ(5) * t229 - pkin(3);
t192 = -t233 * t378 + t298;
t214 = t377 * t229;
t140 = t192 * t232 + t214 * t228;
t433 = -qJD(6) * t140 - t228 * t438 + t232 * t437;
t139 = -t192 * t228 + t214 * t232;
t432 = qJD(6) * t139 + t228 * t437 + t232 * t438;
t308 = t229 * t339;
t288 = qJD(4) * t308;
t300 = qJD(3) * t338;
t143 = -t226 * t317 - t233 * t300 + t288;
t319 = qJD(3) * t230;
t301 = t224 * t319;
t265 = t378 * t301;
t431 = -pkin(5) * t143 - t265 + t435;
t189 = t226 * t229 + t233 * t339;
t289 = t229 * t300;
t144 = qJD(4) * t189 + t289;
t329 = t231 * t234;
t330 = t230 * t235;
t261 = t226 * t329 + t330;
t149 = t261 * t323;
t182 = t191 * qJD(3);
t246 = qJ(5) * t143 - qJD(5) * t189 + t182;
t430 = -t144 * t378 + t149 - t246;
t332 = t229 * t234;
t83 = t233 * t114 + t229 * t179;
t429 = -(-pkin(5) * t332 + qJ(5) * t230) * t321 - t83 - t377 * t318;
t299 = qJD(3) * t321;
t286 = t234 * t299;
t131 = qJD(2) * t288 - qJD(4) * t264 - t233 * t286;
t287 = t230 * t299;
t108 = -t131 * mrSges(6,1) + mrSges(6,2) * t287;
t138 = (t258 + t305) * t321;
t304 = t226 * t322;
t247 = (qJD(2) * t283 - t202) * t224 + t304;
t252 = t259 * qJD(2);
t290 = t227 * t300;
t71 = (t202 * t335 - t184) * qJD(3) + (t225 * t252 + t290) * qJD(1);
t426 = qJD(4) * t247 + t71;
t99 = pkin(10) * t296 + t115;
t17 = t138 * t233 - t229 * t426 - t99 * t317;
t12 = -pkin(4) * t287 - t17;
t428 = m(6) * t12 + t108;
t98 = -pkin(3) * t296 - t114;
t244 = -t165 * qJ(5) + t98;
t47 = t164 * pkin(4) + t244;
t45 = t229 * t99 - t233 * t247;
t263 = pkin(5) * t165 + t45;
t422 = qJD(5) + t263;
t27 = -t378 * t441 + t422;
t38 = t164 * t378 + t244;
t7 = -t228 * t38 + t232 * t27;
t8 = t228 * t27 + t232 * t38;
t427 = t7 * mrSges(7,1) + t98 * mrSges(5,2) - t8 * mrSges(7,2) - t47 * mrSges(6,3) + Ifges(5,5) * t362 + t440;
t355 = -qJD(3) / 0.2e1;
t113 = -mrSges(6,2) * t164 - mrSges(6,3) * t165;
t326 = -mrSges(4,1) * t296 + mrSges(5,1) * t164 + mrSges(5,2) * t165 + mrSges(4,3) * t303;
t425 = -t113 - t326;
t109 = mrSges(5,1) * t287 + mrSges(5,3) * t131;
t424 = -m(5) * t17 - t109 + t428;
t46 = t229 * t247 + t233 * t99;
t40 = -qJ(5) * t441 - t46;
t423 = t40 * mrSges(6,1) - t46 * mrSges(5,3);
t402 = -qJD(5) - t45;
t39 = -pkin(4) * t441 - t402;
t157 = Ifges(6,6) * t164;
t93 = Ifges(6,4) * t441 - t165 * Ifges(6,2) + t157;
t421 = t93 / 0.2e1 - t45 * mrSges(5,3) - t39 * mrSges(6,1);
t420 = -t98 * mrSges(5,1) + t47 * mrSges(6,2) + t441 * t442 + Ifges(5,6) * t362 + (Ifges(6,3) + Ifges(5,2)) * t367 + (Ifges(6,6) + Ifges(5,4)) * t364;
t365 = -t165 / 0.2e1;
t366 = t164 / 0.2e1;
t406 = Ifges(6,4) - Ifges(5,5);
t419 = -Ifges(6,2) * t365 - Ifges(6,6) * t366 - t406 * t362 + t427 + t440;
t418 = -Ifges(6,4) / 0.2e1;
t10 = -pkin(5) * t131 - qJD(2) * t265 - t17;
t132 = qJD(4) * t165 + t229 * t286;
t253 = t261 * qJD(2);
t291 = t227 * t301;
t72 = t399 * qJD(3) + (t225 * t253 + t291) * qJD(1);
t238 = qJ(5) * t131 - qJD(5) * t165 + t72;
t18 = t132 * t378 + t238;
t1 = qJD(6) * t7 + t10 * t228 + t18 * t232;
t417 = t1 * mrSges(7,2);
t2 = -qJD(6) * t8 + t10 * t232 - t18 * t228;
t416 = t2 * mrSges(7,1);
t371 = -t131 / 0.2e1;
t415 = -t132 / 0.2e1;
t414 = t303 / 0.2e1;
t413 = Ifges(4,6) * t355;
t116 = -t229 * t174 + t175 * t233;
t104 = pkin(4) * t338 - t116;
t70 = pkin(5) * t189 + pkin(11) * t338 + t104;
t188 = -t233 * t226 + t308;
t173 = t218 + (-pkin(2) * t234 - pkin(3)) * t226;
t251 = -qJ(5) * t189 + t173;
t81 = t188 * t378 + t251;
t25 = -t228 * t81 + t232 * t70;
t412 = qJD(6) * t25 + t228 * t431 - t232 * t430;
t26 = t228 * t70 + t232 * t81;
t411 = -qJD(6) * t26 + t228 * t430 + t232 * t431;
t405 = Ifges(6,5) - Ifges(5,6);
t134 = -mrSges(5,2) * t441 - mrSges(5,3) * t164;
t136 = mrSges(6,1) * t164 - mrSges(6,3) * t441;
t325 = t134 - t136;
t64 = -mrSges(7,1) * t124 + mrSges(7,2) * t125;
t403 = -t64 - t325;
t135 = mrSges(5,1) * t441 - mrSges(5,3) * t165;
t137 = mrSges(6,1) * t165 + mrSges(6,2) * t441;
t324 = t135 - t137;
t401 = t149 - t182;
t400 = -t150 + t181;
t398 = t226 * t327 - t331;
t397 = t1 * t228 + t2 * t232;
t312 = Ifges(5,5) / 0.2e1 + t418;
t396 = t312 * t441 - t421 + t427 + t434;
t255 = t174 * t318 - t175 * t317 - t229 * t180 - t233 * t181;
t41 = -t224 * (qJ(5) * t319 - qJD(5) * t234) + t255;
t107 = mrSges(6,1) * t132 - mrSges(6,3) * t287;
t16 = t229 * t138 + t233 * t426 - t318 * t99;
t11 = -qJ(5) * t287 - qJD(5) * t441 - t16;
t110 = -mrSges(5,2) * t287 - mrSges(5,3) * t132;
t393 = m(5) * t16 - m(6) * t11 - t107 + t110;
t156 = -t202 * t224 + t304;
t310 = Ifges(5,6) / 0.2e1 + t442;
t313 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t349 = Ifges(4,6) * t226;
t392 = -t310 * t164 + t312 * t165 + t313 * t441 + t156 * mrSges(4,1) + t39 * mrSges(6,2) + t413 - (t349 + (Ifges(4,4) * t230 + Ifges(4,2) * t234) * t224) * qJD(2) / 0.2e1 + Ifges(5,5) * t364 + Ifges(5,6) * t367 + Ifges(6,4) * t365 + Ifges(6,5) * t366 - t115 * mrSges(4,3) - t40 * mrSges(6,3) - t45 * mrSges(5,1) - t46 * mrSges(5,2) + (Ifges(5,3) + Ifges(6,1)) * t362;
t352 = Ifges(7,4) * t228;
t273 = Ifges(7,2) * t232 + t352;
t351 = Ifges(7,4) * t232;
t275 = Ifges(7,1) * t228 + t351;
t278 = mrSges(7,1) * t232 - mrSges(7,2) * t228;
t280 = t228 * t7 - t232 * t8;
t357 = t164 * pkin(5);
t30 = -t40 - t357;
t348 = Ifges(7,6) * t232;
t350 = Ifges(7,5) * t228;
t360 = -t232 / 0.2e1;
t361 = -t228 / 0.2e1;
t369 = -t159 / 0.2e1;
t373 = -t125 / 0.2e1;
t375 = -t124 / 0.2e1;
t353 = Ifges(7,4) * t125;
t43 = Ifges(7,2) * t124 + Ifges(7,6) * t159 + t353;
t123 = Ifges(7,4) * t124;
t44 = Ifges(7,1) * t125 + Ifges(7,5) * t159 + t123;
t390 = mrSges(7,3) * t280 + (t348 + t350) * t369 + t273 * t375 + t275 * t373 + t30 * t278 + t360 * t43 + t361 * t44;
t389 = -Ifges(5,4) * t364 - Ifges(5,2) * t367 + Ifges(6,6) * t365 + Ifges(6,3) * t366 + t362 * t405 - t420;
t55 = qJD(6) * t124 + t132 * t228 + t232 * t287;
t56 = -qJD(6) * t125 + t132 * t232 - t228 * t287;
t15 = t55 * Ifges(7,1) + t56 * Ifges(7,4) - t131 * Ifges(7,5);
t387 = t15 / 0.2e1;
t386 = t43 / 0.2e1;
t385 = t44 / 0.2e1;
t384 = t55 / 0.2e1;
t383 = t56 / 0.2e1;
t382 = Ifges(6,2) * t371 + Ifges(6,6) * t415 + t287 * t418;
t52 = Ifges(7,5) * t55;
t51 = Ifges(7,6) * t56;
t354 = qJD(3) / 0.2e1;
t141 = -t225 * t398 - t227 * t338;
t344 = t72 * t141;
t343 = -t136 + t64;
t341 = qJ(5) * t164;
t334 = t228 * t229;
t333 = t229 * t232;
t328 = t232 * t234;
t117 = t233 * t174 + t229 * t175;
t13 = -Ifges(7,3) * t131 + t51 + t52;
t314 = -Ifges(5,1) / 0.2e1 - Ifges(6,2) / 0.2e1;
t311 = -Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t309 = Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1;
t216 = Ifges(4,4) * t302;
t295 = t156 * mrSges(4,2) + Ifges(4,1) * t414 + Ifges(4,5) * t296 / 0.2e1 + t216 / 0.2e1;
t292 = t225 * t231 * t321;
t285 = t416 - t417;
t284 = -t72 * mrSges(4,1) - t71 * mrSges(4,2);
t279 = -mrSges(4,1) * t234 + mrSges(4,2) * t230;
t277 = mrSges(7,1) * t228 + mrSges(7,2) * t232;
t276 = Ifges(7,1) * t232 - t352;
t274 = -Ifges(7,2) * t228 + t351;
t272 = Ifges(7,5) * t232 - Ifges(7,6) * t228;
t33 = -mrSges(7,1) * t131 - mrSges(7,3) * t55;
t34 = mrSges(7,2) * t131 + mrSges(7,3) * t56;
t269 = t228 * t34 + t232 * t33;
t87 = -mrSges(7,2) * t159 + mrSges(7,3) * t124;
t88 = mrSges(7,1) * t159 - mrSges(7,3) * t125;
t268 = -t228 * t88 + t232 * t87;
t267 = t229 * t40 + t233 * t39;
t266 = -t229 * t46 + t233 * t45;
t260 = t226 * t330 + t329;
t142 = t225 * t260 + t227 * t339;
t187 = -t224 * t225 * t235 + t226 * t227;
t101 = t142 * t229 - t233 * t187;
t53 = t101 * t232 - t141 * t228;
t54 = t101 * t228 + t141 * t232;
t102 = t142 * t233 + t187 * t229;
t103 = qJ(5) * t338 - t117;
t262 = -t188 * t228 + t224 * t328;
t145 = t188 * t232 + t228 * t338;
t254 = t398 * qJD(3);
t250 = -t17 * mrSges(5,1) + t16 * mrSges(5,2) - t12 * mrSges(6,2) + t11 * mrSges(6,3);
t243 = (t254 + t252) * t225;
t241 = t309 * t165 + t310 * t441 + t420 - t423;
t240 = t243 + t290;
t211 = Ifges(6,1) * t287;
t210 = Ifges(4,5) * t286;
t209 = Ifges(5,3) * t287;
t206 = -pkin(4) * t233 + t298;
t186 = -qJ(5) * t317 + t297;
t178 = t279 * t321;
t177 = -mrSges(4,2) * t296 + mrSges(4,3) * t302;
t171 = (mrSges(4,1) * t230 + mrSges(4,2) * t234) * t299;
t155 = (t228 * t332 + t230 * t232) * t321;
t154 = (-t228 * t230 + t229 * t328) * t321;
t130 = Ifges(6,4) * t131;
t129 = Ifges(5,5) * t131;
t128 = Ifges(6,5) * t132;
t127 = Ifges(5,6) * t132;
t122 = t150 * t233 + t229 * t293;
t111 = pkin(4) * t165 + t341;
t100 = pkin(4) * t188 + t251;
t97 = t291 + (qJD(3) * t260 + t253) * t225;
t86 = -t302 * t340 + t294;
t85 = t165 * t378 + t341;
t84 = -pkin(5) * t188 - t103;
t78 = qJD(6) * t262 + t144 * t232 - t228 * t301;
t77 = qJD(6) * t145 + t144 * t228 + t232 * t301;
t74 = mrSges(5,1) * t132 - mrSges(5,2) * t131;
t73 = -mrSges(6,2) * t132 + mrSges(6,3) * t131;
t69 = -pkin(4) * t303 - t82;
t65 = -qJ(5) * t303 - t83;
t63 = -t131 * Ifges(5,1) - t132 * Ifges(5,4) + Ifges(5,5) * t287;
t62 = -t131 * Ifges(5,4) - t132 * Ifges(5,2) + Ifges(5,6) * t287;
t60 = Ifges(6,5) * t287 + t131 * Ifges(6,6) + t132 * Ifges(6,3);
t50 = pkin(4) * t144 + t246;
t48 = -pkin(4) * t301 - t59;
t37 = t46 - t357;
t32 = qJD(4) * t102 + t227 * t289 + t229 * t243 - t233 * t292;
t29 = -pkin(5) * t144 - t41;
t24 = pkin(4) * t132 + t238;
t21 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t20 = t228 * t37 + t232 * t85;
t19 = -t228 * t85 + t232 * t37;
t14 = t55 * Ifges(7,4) + t56 * Ifges(7,2) - t131 * Ifges(7,6);
t9 = -pkin(5) * t132 - t11;
t6 = -qJD(6) * t54 - t228 * t97 + t232 * t32;
t5 = qJD(6) * t53 + t228 * t32 + t232 * t97;
t3 = [m(7) * (t1 * t54 + t2 * t53 + t5 * t8 + t6 * t7) + m(5) * t344 + t240 * t177 + m(4) * (t115 * t290 + t344 + t71 * t142 + (t115 * t254 + (t115 * t259 + (qJD(1) * t187 + t156) * t231 * t224) * qJD(2)) * t225) - t142 * mrSges(4,3) * t287 + t53 * t33 + t54 * t34 + t5 * t87 + t6 * t88 + t178 * t292 + t187 * t171 + (-mrSges(3,1) * t231 - mrSges(3,2) * t235) * t225 * qJD(2) ^ 2 + (m(5) * t45 + m(6) * t39 - t324) * t32 + t424 * t101 + (-m(4) * t114 + m(5) * t98 + m(6) * t47 - t425) * t97 + (-m(5) * t46 + m(6) * t40 - m(7) * t30 + t403) * (t142 * t318 - t187 * t317 - t229 * t292 - t233 * t240) + (m(6) * t24 + mrSges(4,3) * t286 + t73 + t74) * t141 + (m(7) * t9 + t21 + t393) * t102; (Ifges(7,5) * t77 + Ifges(7,6) * t78) * t368 - t255 * t134 - t262 * t387 + (t1 * t145 + t2 * t262 - t7 * t77 + t78 * t8) * mrSges(7,3) + (-Ifges(7,4) * t262 + Ifges(7,2) * t145 + Ifges(7,6) * t189) * t383 + (-Ifges(7,1) * t262 + Ifges(7,4) * t145 + Ifges(7,5) * t189) * t384 + t9 * (-mrSges(7,1) * t145 - mrSges(7,2) * t262) + (-Ifges(5,4) * t188 - Ifges(7,5) * t262 + Ifges(7,6) * t145 + (Ifges(5,1) + Ifges(7,3)) * t189) * t371 + (t11 * t188 + t12 * t189) * mrSges(6,1) + (-t16 * t188 - t17 * t189) * mrSges(5,3) + t425 * t149 + (Ifges(7,4) * t77 + Ifges(7,2) * t78) * t374 + t131 * (-Ifges(6,2) * t189 + Ifges(6,6) * t188) / 0.2e1 + (t389 + t423) * t144 + t189 * t382 + t77 * t385 + t78 * t386 + (Ifges(7,1) * t77 + Ifges(7,4) * t78) * t372 + (t63 + t13) * t189 / 0.2e1 + (t116 * t17 + t117 * t16 + t173 * t72 - t401 * t98 + (-t122 - t255) * t46 + t435 * t45) * m(5) + (-t419 + t421) * t143 + t324 * t121 + t326 * t182 + (t100 * t24 + t103 * t11 + t104 * t12 + (-t149 + t50) * t47 + (t122 + t41) * t40 + (-t121 + t48) * t39) * m(6) + t132 * (-Ifges(6,6) * t189 + Ifges(6,3) * t188) / 0.2e1 + t25 * t33 + t26 * t34 + t403 * t122 + t400 * t177 + (t114 * t401 + t115 * t400 - t190 * t72 + t191 * t71) * m(4) + t29 * t64 + t30 * (-mrSges(7,1) * t78 + mrSges(7,2) * t77) + t84 * t21 + t411 * t88 + t412 * t87 + (t1 * t26 + t2 * t25 + t84 * t9 + t412 * t8 + t411 * t7 + (-t122 + t29) * t30) * m(7) + (t72 * mrSges(4,3) * t230 - pkin(2) * t171 + (t71 * mrSges(4,3) + t129 / 0.2e1 + t127 / 0.2e1 - t209 / 0.2e1 - t211 / 0.2e1 - t130 / 0.2e1 - t128 / 0.2e1 + t310 * t132 + t312 * t131 + t250) * t234 + (-m(4) * t156 - t178 + (-m(4) * pkin(2) + t279) * t321) * t305 + (((t354 + t320) * Ifges(4,5) + (-qJD(2) * t190 - t114) * mrSges(4,3) + t295) * t234 + (t413 + (-0.3e1 / 0.2e1 * t349 - t191 * mrSges(4,3) + t312 * t189 - t310 * t188 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - t313) * t338) * qJD(2) + t392) * t230) * qJD(3) + (-0.3e1 / 0.2e1 * t230 ^ 2 + 0.3e1 / 0.2e1 * t234 ^ 2) * Ifges(4,4) * t299) * t224 + t100 * t73 - t189 * t417 + t103 * t107 + t104 * t108 + t50 * t113 + t116 * t109 + t117 * t110 + t59 * t135 + t41 * t136 + t48 * t137 + t145 * t14 / 0.2e1 + t173 * t74 + t188 * t60 / 0.2e1 - t188 * t62 / 0.2e1 + t72 * (mrSges(5,1) * t188 + mrSges(5,2) * t189) + t24 * (-mrSges(6,2) * t188 - mrSges(6,3) * t189) + (t210 / 0.2e1 + t284) * t226 + (Ifges(5,4) * t189 - Ifges(5,2) * t188) * t415 + t189 * t416; (t62 / 0.2e1 - t60 / 0.2e1 - t72 * mrSges(5,1) + t24 * mrSges(6,2) - t11 * mrSges(6,1) + t16 * mrSges(5,3) + t15 * t361 + t14 * t360 + t9 * t278 - t55 * t275 / 0.2e1 - t56 * t273 / 0.2e1 + t311 * t132 + (t350 / 0.2e1 + t348 / 0.2e1 - t309) * t131 + (-t1 * t232 + t2 * t228) * mrSges(7,3) + t393 * pkin(10) + (t44 * t360 + t228 * t386 + t274 * t375 + t276 * t373 + t272 * t369 - t30 * t277 + (t228 * t8 + t232 * t7) * mrSges(7,3)) * qJD(6)) * t233 + (-t86 + t186) * t113 + (-t8 * t154 + t7 * t155) * mrSges(7,3) + t284 + m(6) * (t186 * t47 + t206 * t24) + t432 * t87 + t433 * t88 + (t1 * t140 + t139 * t2 + t215 * t9 + t30 * t429 + t432 * t8 + t433 * t7) * m(7) + (t30 * (-mrSges(7,1) * t333 + mrSges(7,2) * t334) + (Ifges(7,4) * t334 + Ifges(7,2) * t333) * t374 + (Ifges(7,1) * t334 + Ifges(7,4) * t333) * t372 + (Ifges(7,5) * t334 + Ifges(7,6) * t333) * t368 + t333 * t386 + t334 * t385 + t266 * mrSges(5,3) + t267 * mrSges(6,1) + (m(5) * t266 + m(6) * t267) * pkin(10) + (t333 * t8 - t334 * t7) * mrSges(7,3) + (-t325 * pkin(10) + t389) * t229 + (-t93 / 0.2e1 - t324 * pkin(10) + t419) * t233) * qJD(4) + t429 * t64 - m(5) * (t115 * t98 - t45 * t82 + t46 * t83) - m(6) * (t39 * t69 + t40 * t65 + t47 * t86) + (-t17 * mrSges(5,3) + t12 * mrSges(6,1) + t13 / 0.2e1 + t382 + t63 / 0.2e1 - t24 * mrSges(6,3) + t72 * mrSges(5,2) + t52 / 0.2e1 + t51 / 0.2e1 - t309 * t132 + (-Ifges(7,3) / 0.2e1 + t314) * t131 + t424 * pkin(10) + t285) * t229 + t210 + (Ifges(7,4) * t155 + Ifges(7,2) * t154) * t375 + (Ifges(7,5) * t155 + Ifges(7,6) * t154) * t369 + (Ifges(7,1) * t155 + Ifges(7,4) * t154) * t373 - t326 * t115 + (-m(5) * t72 - t74) * pkin(3) + (((t355 + t320 / 0.2e1) * Ifges(4,6) + Ifges(4,4) * t414 + (-t229 * t406 - t233 * t405) * t354 - t392) * t230 + (-t216 / 0.2e1 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t339 - t226 * Ifges(4,5) / 0.2e1) * qJD(2) + Ifges(4,5) * t355 + t114 * mrSges(4,3) + (t164 * t311 + t241) * t229 + (t164 * t309 + t165 * t314 - t396) * t233 - t295) * t234) * t321 - t83 * t134 - t82 * t135 - t65 * t136 - t69 * t137 + t139 * t33 + t140 * t34 - t154 * t43 / 0.2e1 - t155 * t44 / 0.2e1 - t30 * (-mrSges(7,1) * t154 + mrSges(7,2) * t155) - t114 * t177 + t206 * t73 + t215 * t21; (t21 - t107) * qJ(5) + t130 + t128 - t129 + (t241 + (t311 - t314) * t164 + t390) * t165 + t390 * qJD(6) + t263 * t64 - t397 * mrSges(7,3) - ((-m(7) * t280 + t268) * qJD(6) + m(7) * t397 + t269) * t378 - t127 + (qJ(5) * t9 - t19 * t7 - t20 * t8 + t30 * t422) * m(7) - t250 + t209 + t211 + t274 * t383 + t276 * t384 + t232 * t387 + t14 * t361 + t272 * t371 + t343 * qJD(5) + t324 * t46 + t325 * t45 + (-t157 / 0.2e1 + t436 + t396) * t164 + (-pkin(4) * t12 - qJ(5) * t11 - t111 * t47 - t39 * t46 + t40 * t402) * m(6) - t20 * t87 - t19 * t88 - pkin(4) * t108 - t111 * t113 + t9 * t277; -t343 * t441 + t268 * qJD(6) + (t113 + t268) * t165 - m(6) * (-t165 * t47 - t40 * t441) + t269 + (-t159 * t280 - t30 * t441 + t397) * m(7) + t428; -t30 * (mrSges(7,1) * t125 + mrSges(7,2) * t124) + (Ifges(7,1) * t124 - t353) * t373 + t43 * t372 + (Ifges(7,5) * t124 - Ifges(7,6) * t125) * t369 - t7 * t87 + t8 * t88 + (t124 * t7 + t125 * t8) * mrSges(7,3) + t285 + t13 + (-Ifges(7,2) * t125 + t123 + t44) * t375;];
tauc  = t3(:);

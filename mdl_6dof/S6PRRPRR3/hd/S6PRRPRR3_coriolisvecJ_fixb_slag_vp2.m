% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:20
% EndTime: 2019-03-08 22:03:02
% DurationCPUTime: 19.86s
% Computational Cost: add. (12272->709), mult. (36844->1047), div. (0->0), fcn. (30220->14), ass. (0->334)
t235 = cos(pkin(7));
t243 = cos(qJ(3));
t244 = cos(qJ(2));
t331 = t243 * t244;
t239 = sin(qJ(3));
t240 = sin(qJ(2));
t336 = t239 * t240;
t259 = -t235 * t336 + t331;
t234 = sin(pkin(6));
t328 = qJD(1) * t234;
t182 = t259 * t328;
t322 = qJD(3) * t243;
t307 = t235 * t322;
t223 = pkin(2) * t307;
t233 = sin(pkin(7));
t370 = -pkin(9) - qJ(4);
t303 = t370 * t239;
t320 = qJD(4) * t243;
t436 = -t223 - (qJD(3) * t303 + t320) * t233 + t182;
t334 = t240 * t243;
t335 = t239 * t244;
t261 = -t235 * t334 - t335;
t181 = t261 * t328;
t341 = t235 * t239;
t227 = pkin(2) * t341;
t321 = qJD(4) * t239;
t342 = t233 * t243;
t435 = -t233 * t321 + (t342 * t370 - t227) * qJD(3) - t181;
t311 = t240 * t328;
t290 = t233 * t311;
t323 = qJD(3) * t239;
t308 = t233 * t323;
t434 = pkin(3) * t308 - t290;
t232 = sin(pkin(13));
t347 = cos(pkin(13));
t432 = -t435 * t232 + t347 * t436;
t296 = t347 * t239;
t195 = (t232 * t243 + t296) * t233;
t190 = qJD(3) * t195;
t295 = t347 * t243;
t255 = -t232 * t239 + t295;
t250 = t233 * t255;
t191 = qJD(3) * t250;
t433 = pkin(4) * t190 - pkin(10) * t191 + t434;
t188 = qJD(2) * t250;
t409 = -t188 + qJD(5);
t340 = t235 * t243;
t228 = pkin(2) * t340;
t170 = pkin(3) * t235 + t233 * t303 + t228;
t205 = pkin(9) * t342 + t227;
t187 = qJ(4) * t342 + t205;
t132 = t232 * t170 + t347 * t187;
t117 = pkin(10) * t235 + t132;
t343 = t233 * t239;
t194 = t232 * t343 - t233 * t295;
t210 = (-pkin(3) * t243 - pkin(2)) * t233;
t144 = pkin(4) * t194 - pkin(10) * t195 + t210;
t238 = sin(qJ(5));
t242 = cos(qJ(5));
t318 = qJD(5) * t242;
t319 = qJD(5) * t238;
t419 = -t117 * t319 + t144 * t318 + t433 * t238 - t242 * t432;
t416 = t232 * t436 + t435 * t347;
t378 = pkin(3) * t232;
t230 = pkin(10) + t378;
t324 = qJD(2) * t243;
t309 = t233 * t324;
t326 = qJD(2) * t233;
t189 = -t232 * t309 - t296 * t326;
t310 = t239 * t326;
t294 = pkin(3) * t310;
t138 = -pkin(4) * t189 - pkin(10) * t188 + t294;
t215 = qJD(2) * pkin(2) + t244 * t328;
t207 = pkin(9) * t326 + t311;
t332 = t243 * t207;
t263 = -t215 * t341 - t332;
t236 = cos(pkin(6));
t327 = qJD(1) * t236;
t143 = (qJ(4) * t324 + t239 * t327) * t233 - t263;
t123 = t232 * t143;
t312 = t233 * t327;
t217 = t243 * t312;
t329 = t215 * t340 + t217;
t142 = (-qJ(4) * t326 - t207) * t239 + t329;
t76 = t142 * t347 - t123;
t44 = t238 * t138 + t242 * t76;
t431 = -t230 * t319 - t44;
t430 = -pkin(11) * t190 - t419;
t167 = t195 * t238 - t242 * t235;
t129 = -qJD(5) * t167 + t191 * t242;
t168 = t195 * t242 + t235 * t238;
t130 = qJD(5) * t168 + t191 * t238;
t429 = pkin(5) * t130 - pkin(11) * t129 - t416;
t297 = t347 * t143;
t75 = t142 * t232 + t297;
t428 = -t75 + t409 * (pkin(5) * t238 - pkin(11) * t242);
t427 = -pkin(11) * t189 - t431;
t414 = t242 * t117 + t238 * t144;
t420 = -qJD(5) * t414 + t432 * t238 + t242 * t433;
t241 = cos(qJ(6));
t237 = sin(qJ(6));
t338 = t237 * t242;
t145 = -t188 * t338 - t189 * t241;
t426 = t237 * t318 + t145;
t225 = qJD(2) * t235 + qJD(3);
t154 = t189 * t238 + t225 * t242;
t180 = qJD(2) * t191;
t108 = qJD(5) * t154 + t180 * t242;
t155 = -t189 * t242 + t225 * t238;
t185 = -t255 * t326 + qJD(5);
t110 = -t155 * t237 + t185 * t241;
t179 = qJD(2) * t190;
t53 = qJD(6) * t110 + t108 * t241 + t179 * t237;
t399 = t53 / 0.2e1;
t111 = t155 * t241 + t185 * t237;
t54 = -qJD(6) * t111 - t108 * t237 + t179 * t241;
t398 = t54 / 0.2e1;
t113 = pkin(3) * t225 + t142;
t68 = t232 * t113 + t297;
t66 = pkin(10) * t225 + t68;
t224 = t235 * t327;
t159 = qJD(4) + t224 + (-pkin(3) * t324 - t215) * t233;
t89 = -pkin(4) * t188 + pkin(10) * t189 + t159;
t35 = t238 * t89 + t242 * t66;
t29 = pkin(11) * t185 + t35;
t67 = t113 * t347 - t123;
t65 = -t225 * pkin(4) - t67;
t36 = -t154 * pkin(5) - t155 * pkin(11) + t65;
t11 = -t237 * t29 + t241 * t36;
t109 = qJD(5) * t155 + t180 * t238;
t150 = t239 * (t215 * t235 + t312) + t332;
t411 = -t150 * qJD(3) + ((-qJ(4) * t322 - t321) * t233 + t181) * qJD(2);
t325 = qJD(2) * t234;
t302 = qJD(1) * t325;
t314 = qJD(3) * t217 + t215 * t307 + t302 * t331;
t83 = -t207 * t323 + (-t311 * t341 + (-qJ(4) * t323 + t320) * t233) * qJD(2) + t314;
t40 = t232 * t83 - t347 * t411;
t19 = pkin(5) * t109 - pkin(11) * t108 + t40;
t41 = t232 * t411 + t347 * t83;
t301 = qJD(3) * t326;
t286 = t239 * t301;
t287 = t240 * t302;
t186 = pkin(3) * t286 + t233 * t287;
t98 = pkin(4) * t179 - pkin(10) * t180 + t186;
t9 = t238 * t98 + t242 * t41 + t89 * t318 - t319 * t66;
t7 = pkin(11) * t179 + t9;
t1 = qJD(6) * t11 + t19 * t237 + t241 * t7;
t425 = t1 * mrSges(7,2);
t12 = t237 * t36 + t241 * t29;
t2 = -qJD(6) * t12 + t19 * t241 - t237 * t7;
t424 = t2 * mrSges(7,1);
t392 = t109 / 0.2e1;
t60 = pkin(11) * t194 + t414;
t131 = t170 * t347 - t232 * t187;
t116 = -t235 * pkin(4) - t131;
t71 = t167 * pkin(5) - t168 * pkin(11) + t116;
t25 = t237 * t71 + t241 * t60;
t423 = -qJD(6) * t25 + t237 * t430 + t429 * t241;
t24 = -t237 * t60 + t241 * t71;
t422 = qJD(6) * t24 + t429 * t237 - t241 * t430;
t18 = -mrSges(7,1) * t54 + mrSges(7,2) * t53;
t84 = mrSges(6,1) * t179 - mrSges(6,3) * t108;
t369 = t18 - t84;
t421 = -pkin(5) * t190 - t420;
t313 = t347 * pkin(3);
t231 = -t313 - pkin(4);
t206 = -t242 * pkin(5) - t238 * pkin(11) + t231;
t165 = t206 * t241 - t230 * t338;
t418 = qJD(6) * t165 + t237 * t428 - t241 * t427;
t333 = t241 * t242;
t166 = t206 * t237 + t230 * t333;
t417 = -qJD(6) * t166 + t237 * t427 + t241 * t428;
t367 = mrSges(5,3) * t189;
t348 = mrSges(5,1) * t225 + mrSges(6,1) * t154 - mrSges(6,2) * t155 + t367;
t317 = qJD(6) * t238;
t413 = t241 * t317 + t426;
t146 = t188 * t333 - t189 * t237;
t412 = t237 * t317 - t241 * t318 + t146;
t10 = -qJD(5) * t35 - t238 * t41 + t242 * t98;
t410 = -t10 * t238 + t242 * t9;
t283 = t1 * t241 - t2 * t237;
t270 = t11 * t241 + t12 * t237;
t272 = Ifges(7,5) * t241 - Ifges(7,6) * t237;
t359 = Ifges(7,4) * t241;
t275 = -Ifges(7,2) * t237 + t359;
t360 = Ifges(7,4) * t237;
t278 = Ifges(7,1) * t241 - t360;
t34 = -t238 * t66 + t242 * t89;
t28 = -pkin(5) * t185 - t34;
t280 = mrSges(7,1) * t237 + mrSges(7,2) * t241;
t106 = Ifges(7,4) * t110;
t153 = qJD(6) - t154;
t47 = Ifges(7,1) * t111 + Ifges(7,5) * t153 + t106;
t350 = t241 * t47;
t386 = t153 / 0.2e1;
t388 = t111 / 0.2e1;
t390 = t110 / 0.2e1;
t361 = Ifges(7,4) * t111;
t46 = Ifges(7,2) * t110 + Ifges(7,6) * t153 + t361;
t400 = -t46 / 0.2e1;
t408 = -t270 * mrSges(7,3) + t28 * t280 + t272 * t386 + t275 * t390 + t278 * t388 + t350 / 0.2e1 + t237 * t400;
t407 = t11 * mrSges(7,1) - t12 * mrSges(7,2);
t353 = t153 * Ifges(7,3);
t354 = t111 * Ifges(7,5);
t355 = t110 * Ifges(7,6);
t45 = t353 + t354 + t355;
t351 = t185 * Ifges(6,6);
t364 = Ifges(6,4) * t155;
t80 = t154 * Ifges(6,2) + t351 + t364;
t406 = t80 / 0.2e1 - t45 / 0.2e1 - t407;
t405 = qJD(2) ^ 2;
t14 = t53 * Ifges(7,4) + t54 * Ifges(7,2) + t109 * Ifges(7,6);
t404 = t14 / 0.2e1;
t403 = Ifges(7,1) * t399 + Ifges(7,4) * t398 + Ifges(7,5) * t392;
t401 = t45 / 0.2e1;
t397 = -t80 / 0.2e1;
t151 = Ifges(6,4) * t154;
t352 = t185 * Ifges(6,5);
t81 = t155 * Ifges(6,1) + t151 + t352;
t395 = -t81 / 0.2e1;
t394 = t108 / 0.2e1;
t393 = -t109 / 0.2e1;
t391 = -t110 / 0.2e1;
t389 = -t111 / 0.2e1;
t387 = -t153 / 0.2e1;
t385 = -t167 / 0.2e1;
t384 = t168 / 0.2e1;
t382 = -t189 / 0.2e1;
t380 = t235 / 0.2e1;
t379 = -t239 / 0.2e1;
t8 = -pkin(5) * t179 - t10;
t375 = t238 * t8;
t373 = t34 * mrSges(6,3);
t372 = t35 * mrSges(6,3);
t368 = mrSges(5,3) * t188;
t366 = Ifges(4,4) * t239;
t365 = Ifges(5,4) * t189;
t363 = Ifges(6,4) * t238;
t362 = Ifges(6,4) * t242;
t358 = Ifges(4,5) * t243;
t262 = t235 * t331 - t336;
t316 = t236 * t342;
t162 = t234 * t262 + t316;
t260 = t235 * t335 + t334;
t163 = t234 * t260 + t236 * t343;
t101 = -t162 * t347 + t163 * t232;
t356 = t101 * t40;
t115 = mrSges(6,1) * t185 - mrSges(6,3) * t155;
t62 = -mrSges(7,1) * t110 + mrSges(7,2) * t111;
t349 = t115 - t62;
t346 = t188 * t238;
t345 = t188 * t242;
t339 = t237 * t238;
t337 = t238 * t241;
t148 = -mrSges(5,1) * t188 - mrSges(5,2) * t189;
t282 = -mrSges(4,1) * t243 + mrSges(4,2) * t239;
t330 = t282 * t326 + t148;
t13 = Ifges(7,5) * t53 + Ifges(7,6) * t54 + Ifges(7,3) * t109;
t315 = Ifges(6,5) * t108 - Ifges(6,6) * t109 + Ifges(6,3) * t179;
t305 = t230 * t318;
t147 = t179 * mrSges(5,1) + t180 * mrSges(5,2);
t292 = mrSges(4,3) * t310;
t291 = mrSges(4,3) * t309;
t289 = t233 * t240 * t325;
t288 = t236 * t308;
t281 = mrSges(4,1) * t239 + mrSges(4,2) * t243;
t279 = Ifges(6,1) * t242 - t363;
t277 = Ifges(7,1) * t237 + t359;
t276 = -Ifges(6,2) * t238 + t362;
t274 = Ifges(7,2) * t241 + t360;
t273 = Ifges(6,5) * t242 - Ifges(6,6) * t238;
t271 = Ifges(7,5) * t237 + Ifges(7,6) * t241;
t30 = mrSges(7,1) * t109 - mrSges(7,3) * t53;
t31 = -mrSges(7,2) * t109 + mrSges(7,3) * t54;
t269 = -t237 * t30 + t241 * t31;
t77 = -mrSges(7,2) * t153 + mrSges(7,3) * t110;
t78 = mrSges(7,1) * t153 - mrSges(7,3) * t111;
t268 = -t237 * t77 - t241 * t78;
t102 = t232 * t162 + t163 * t347;
t203 = -t233 * t234 * t244 + t235 * t236;
t87 = t102 * t242 + t203 * t238;
t51 = t101 * t241 - t237 * t87;
t52 = t101 * t237 + t241 * t87;
t43 = t138 * t242 - t238 * t76;
t86 = t102 * t238 - t203 * t242;
t72 = -t117 * t238 + t144 * t242;
t141 = t168 * t241 + t194 * t237;
t140 = -t168 * t237 + t194 * t241;
t256 = t225 * (-Ifges(4,6) * t239 + t358);
t252 = t281 * t326;
t251 = t261 * qJD(2);
t92 = (-qJD(3) * t207 - t235 * t287) * t239 + t314;
t93 = t263 * qJD(3) + (t234 * t251 - t288) * qJD(1);
t249 = t93 * mrSges(4,1) - t40 * mrSges(5,1) - t92 * mrSges(4,2) - t41 * mrSges(5,2);
t222 = Ifges(4,4) * t309;
t219 = t301 * t358;
t204 = -pkin(9) * t343 + t228;
t201 = t205 * qJD(3);
t200 = -pkin(9) * t308 + t223;
t198 = -mrSges(4,2) * t225 + t291;
t197 = mrSges(4,1) * t225 - t292;
t192 = qJD(3) * t252;
t184 = Ifges(5,4) * t188;
t183 = -t215 * t233 + t224;
t176 = Ifges(5,5) * t180;
t175 = Ifges(5,6) * t179;
t172 = Ifges(4,1) * t310 + t225 * Ifges(4,5) + t222;
t171 = t225 * Ifges(4,6) + (t243 * Ifges(4,2) + t366) * t326;
t160 = -mrSges(5,2) * t225 + t368;
t149 = -t207 * t239 + t329;
t135 = -t189 * Ifges(5,1) + t225 * Ifges(5,5) + t184;
t134 = t188 * Ifges(5,2) + t225 * Ifges(5,6) - t365;
t128 = qJD(3) * t316 + (qJD(2) * t259 + qJD(3) * t262) * t234;
t127 = -t288 + (-qJD(3) * t260 + t251) * t234;
t114 = -mrSges(6,2) * t185 + mrSges(6,3) * t154;
t97 = pkin(5) * t155 - pkin(11) * t154;
t85 = -mrSges(6,2) * t179 - mrSges(6,3) * t109;
t79 = t155 * Ifges(6,5) + t154 * Ifges(6,6) + t185 * Ifges(6,3);
t70 = t232 * t127 + t128 * t347;
t69 = -t127 * t347 + t128 * t232;
t61 = mrSges(6,1) * t109 + mrSges(6,2) * t108;
t59 = -pkin(5) * t194 - t72;
t58 = -qJD(6) * t141 - t129 * t237 + t190 * t241;
t57 = qJD(6) * t140 + t129 * t241 + t190 * t237;
t56 = t108 * Ifges(6,1) - t109 * Ifges(6,4) + Ifges(6,5) * t179;
t55 = t108 * Ifges(6,4) - t109 * Ifges(6,2) + Ifges(6,6) * t179;
t37 = pkin(5) * t189 - t43;
t33 = qJD(5) * t87 + t238 * t70 - t242 * t289;
t32 = -qJD(5) * t86 + t238 * t289 + t242 * t70;
t21 = t237 * t97 + t241 * t34;
t20 = -t237 * t34 + t241 * t97;
t6 = -qJD(6) * t52 - t237 * t32 + t241 * t69;
t5 = qJD(6) * t51 + t237 * t69 + t241 * t32;
t3 = [t101 * t61 + t32 * t114 + t127 * t197 + t128 * t198 + t70 * t160 + t51 * t30 + t52 * t31 + t5 * t77 + t6 * t78 + t87 * t85 + t369 * t86 - t348 * t69 - t349 * t33 + (-mrSges(3,1) * t240 - mrSges(3,2) * t244) * t405 * t234 + (t192 + t147) * t203 + (t101 * t180 - t102 * t179) * mrSges(5,3) + (t330 * t240 * t234 + (-t162 * t243 - t163 * t239) * qJD(3) * mrSges(4,3)) * t326 + m(4) * (t127 * t149 + t128 * t150 + t162 * t93 + t163 * t92 + (qJD(1) * t203 + t183) * t289) + m(5) * (t102 * t41 + t159 * t289 + t186 * t203 - t67 * t69 + t68 * t70 + t356) + m(6) * (-t10 * t86 + t32 * t35 - t33 * t34 + t65 * t69 + t87 * t9 + t356) + m(7) * (t1 * t52 + t11 * t6 + t12 * t5 + t2 * t51 + t28 * t33 + t8 * t86); -t432 * t160 + (-t190 * t68 - t191 * t67 - t194 * t41 + t195 * t40) * mrSges(5,3) + m(4) * (-t149 * t201 + t150 * t200 + t204 * t93 + t205 * t92) + (-t201 - t181) * t197 + t414 * t85 + (t10 * t72 + t116 * t40 + t34 * t420 + t35 * t419 + t414 * t9 - t416 * t65) * m(6) + (t219 / 0.2e1 + t176 / 0.2e1 - t175 / 0.2e1 + t249) * t235 + (-Ifges(5,4) * t195 - Ifges(5,6) * t235 / 0.2e1 + Ifges(6,5) * t384 + Ifges(6,6) * t385 - t132 * mrSges(5,3) + (Ifges(5,2) + Ifges(6,3) / 0.2e1) * t194) * t179 + (-pkin(2) * t192 - t330 * t311 + (-t239 * t93 + t243 * t92) * mrSges(4,3) + (t239 * pkin(3) * t148 + t256 / 0.2e1 + t183 * t281 + t243 * t172 / 0.2e1 + t171 * t379 + (-t149 * t243 - t150 * t239) * mrSges(4,3)) * qJD(3)) * t233 + (-t131 * mrSges(5,3) + Ifges(5,1) * t195 - Ifges(5,4) * t194 + Ifges(5,5) * t380) * t180 + ((-m(4) * pkin(2) + t282) * t290 + ((Ifges(4,5) * t380 - t204 * mrSges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,4) * t342) * t243 + (-t205 * mrSges(4,3) - Ifges(4,6) * t235 + (-0.3e1 / 0.2e1 * t366 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t243) * t233) * t239) * qJD(3)) * t326 + t194 * t315 / 0.2e1 + (t200 - t182) * t198 + t225 * (Ifges(5,5) * t191 - Ifges(5,6) * t190) / 0.2e1 + t210 * t147 + t9 * (-mrSges(6,2) * t194 - mrSges(6,3) * t167) + t10 * (mrSges(6,1) * t194 - mrSges(6,3) * t168) + t186 * (mrSges(5,1) * t194 + mrSges(5,2) * t195) - t190 * t134 / 0.2e1 + t188 * (Ifges(5,4) * t191 - Ifges(5,2) * t190) / 0.2e1 + t159 * (mrSges(5,1) * t190 + mrSges(5,2) * t191) + t191 * t135 / 0.2e1 + t185 * (Ifges(6,5) * t129 - Ifges(6,6) * t130 + Ifges(6,3) * t190) / 0.2e1 + t154 * (Ifges(6,4) * t129 - Ifges(6,2) * t130 + Ifges(6,6) * t190) / 0.2e1 + t155 * (Ifges(6,1) * t129 - Ifges(6,4) * t130 + Ifges(6,5) * t190) / 0.2e1 + t35 * (-mrSges(6,2) * t190 - mrSges(6,3) * t130) + t34 * (mrSges(6,1) * t190 - mrSges(6,3) * t129) + t190 * t79 / 0.2e1 + t2 * (mrSges(7,1) * t167 - mrSges(7,3) * t141) + t1 * (-mrSges(7,2) * t167 + mrSges(7,3) * t140) + t167 * t13 / 0.2e1 + t40 * (mrSges(6,1) * t167 + mrSges(6,2) * t168) + t8 * (-mrSges(7,1) * t140 + mrSges(7,2) * t141) + t129 * t81 / 0.2e1 + t12 * (-mrSges(7,2) * t130 + mrSges(7,3) * t58) + t11 * (mrSges(7,1) * t130 - mrSges(7,3) * t57) + t65 * (mrSges(6,1) * t130 + mrSges(6,2) * t129) + t116 * t61 + t72 * t84 + t59 * t18 + t57 * t47 / 0.2e1 + t28 * (-mrSges(7,1) * t58 + mrSges(7,2) * t57) + t58 * t46 / 0.2e1 + (Ifges(7,5) * t141 + Ifges(7,6) * t140 + Ifges(7,3) * t167) * t392 + (Ifges(6,4) * t168 - Ifges(6,2) * t167 + Ifges(6,6) * t194) * t393 + (Ifges(6,1) * t168 - Ifges(6,4) * t167 + Ifges(6,5) * t194) * t394 + t130 * t397 + (Ifges(7,4) * t141 + Ifges(7,2) * t140 + Ifges(7,6) * t167) * t398 + (Ifges(7,1) * t141 + Ifges(7,4) * t140 + Ifges(7,5) * t167) * t399 + t130 * t401 + t141 * t403 + t140 * t404 + t56 * t384 + t55 * t385 + (Ifges(7,5) * t57 + Ifges(7,6) * t58 + Ifges(7,3) * t130) * t386 + (Ifges(7,1) * t57 + Ifges(7,4) * t58 + Ifges(7,5) * t130) * t388 + (Ifges(7,4) * t57 + Ifges(7,2) * t58 + Ifges(7,6) * t130) * t390 + (Ifges(5,1) * t191 - Ifges(5,4) * t190) * t382 + t24 * t30 + t25 * t31 - m(4) * (t149 * t181 + t150 * t182 + t183 * t290) + t416 * t348 + t419 * t114 + t420 * t115 + t421 * t62 + (-t131 * t40 + t132 * t41 + t159 * t434 + t186 * t210 + t416 * t67 - t432 * t68) * m(5) + t422 * t77 + t423 * t78 + (t1 * t25 + t11 * t423 + t12 * t422 + t2 * t24 + t28 * t421 + t59 * t8) * m(7); t431 * t114 + (-t305 - t43) * t115 - t14 * t339 / 0.2e1 + (Ifges(7,1) * t146 + Ifges(7,4) * t145) * t389 + (Ifges(7,5) * t146 + Ifges(7,6) * t145) * t387 + t409 * t65 * (mrSges(6,1) * t238 + mrSges(6,2) * t242) + t34 * mrSges(6,1) * t189 - t35 * mrSges(6,2) * t189 + (t291 - t198) * t149 - t68 * t367 - t154 * (-Ifges(6,6) * t189 + t188 * t276) / 0.2e1 - t155 * (-Ifges(6,5) * t189 + t188 * t279) / 0.2e1 - (t256 + (-Ifges(4,2) * t310 + t172 + t222) * t243) * t326 / 0.2e1 - t318 * t373 + (t292 + t197) * t150 + t249 + t219 + (-t372 + t397 + t401 + t407) * t319 + (Ifges(7,5) * t389 + Ifges(7,6) * t391 + Ifges(7,3) * t387 + t406) * t346 - (t237 * t47 + t241 * t46) * t317 / 0.2e1 + (t81 + t350) * t318 / 0.2e1 + (t154 * t276 + t155 * t279 + t185 * t273) * qJD(5) / 0.2e1 + (Ifges(5,1) * t188 + t365 + t79) * t189 / 0.2e1 - (Ifges(5,2) * t189 + t135 + t184) * t188 / 0.2e1 + (-t37 + t305) * t62 + (-t179 * t378 - t180 * t313) * mrSges(5,3) - t185 * (-Ifges(6,3) * t189 + t188 * t273) / 0.2e1 + (Ifges(7,4) * t146 + Ifges(7,2) * t145) * t391 - t148 * t294 - Ifges(4,6) * t286 + t171 * t310 / 0.2e1 + t233 ^ 2 * t405 * (Ifges(4,1) * t243 - t366) * t379 - t175 + t176 + t179 * (Ifges(6,5) * t238 + Ifges(6,6) * t242) / 0.2e1 + t40 * (-mrSges(6,1) * t242 + mrSges(6,2) * t238) + t242 * t55 / 0.2e1 - t242 * t13 / 0.2e1 + t238 * t56 / 0.2e1 - t225 * (Ifges(5,5) * t188 + Ifges(5,6) * t189) / 0.2e1 + t231 * t61 - t159 * (-mrSges(5,1) * t189 + mrSges(5,2) * t188) + t165 * t30 + t166 * t31 - t76 * t160 - t146 * t47 / 0.2e1 + (-Ifges(7,3) * t242 + t238 * t272) * t392 + (Ifges(6,2) * t242 + t363) * t393 + (Ifges(6,1) * t238 + t362) * t394 + t345 * t395 + (-Ifges(7,6) * t242 + t238 * t275) * t398 + (-Ifges(7,5) * t242 + t238 * t278) * t399 + t337 * t403 + (-t271 * t317 + (Ifges(7,3) * t238 + t242 * t272) * qJD(5)) * t386 + (-t277 * t317 + (Ifges(7,5) * t238 + t242 * t278) * qJD(5)) * t388 + (-t274 * t317 + (Ifges(7,6) * t238 + t242 * t275) * qJD(5)) * t390 + t67 * t368 + t280 * t375 + t134 * t382 - t183 * t252 + t426 * t400 + (t34 * t345 + t346 * t35 + t410) * mrSges(6,3) + (mrSges(7,1) * t413 - mrSges(7,2) * t412) * t28 + (-t1 * t339 + t11 * t412 - t12 * t413 - t2 * t337) * mrSges(7,3) + ((t232 * t41 - t347 * t40) * pkin(3) - t159 * t294 + t67 * t75 - t68 * t76) * m(5) + t348 * t75 + t417 * t78 + t418 * t77 + t242 * t425 - t242 * t424 + (t231 * t40 - t34 * t43 - t35 * t44 - t65 * t75) * m(6) + (t1 * t166 + t417 * t11 + t418 * t12 + t165 * t2 - t28 * t37) * m(7) + (t369 * t238 + t242 * t85 + ((-t238 * t35 - t242 * t34) * qJD(5) + t410) * m(6) + (t28 * t318 + t375) * m(7)) * t230; -t145 * t78 - t146 * t77 - t188 * t160 - t348 * t189 + (-t188 * t114 + (-t237 * t78 + t241 * t77 + t114) * qJD(5) - t369) * t242 + (qJD(6) * t268 - t349 * t409 + t269 + t85) * t238 + t147 + ((-t8 + (-t11 * t237 + t12 * t241) * qJD(5)) * t242 + (qJD(5) * t28 - qJD(6) * t270 + t283) * t238 - t11 * t145 - t12 * t146 - t28 * t346) * m(7) + (t10 * t242 + t189 * t65 + t238 * t9 + t409 * (-t238 * t34 + t242 * t35)) * m(6) + (-t188 * t68 - t189 * t67 + t186) * m(5); t349 * t35 + t277 * t399 + (-t352 / 0.2e1 + t395 - t65 * mrSges(6,2) - t151 / 0.2e1 + t373 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t155 - t408) * t154 + t408 * qJD(6) + t271 * t392 + t274 * t398 + t283 * mrSges(7,3) + t8 * (-mrSges(7,1) * t241 + mrSges(7,2) * t237) + t241 * t404 + t237 * t403 - t34 * t114 - t20 * t78 - t21 * t77 - pkin(5) * t18 - t9 * mrSges(6,2) + t10 * mrSges(6,1) + (t372 - t355 / 0.2e1 + t351 / 0.2e1 - t354 / 0.2e1 - t353 / 0.2e1 - t65 * mrSges(6,1) + t364 / 0.2e1 + t406) * t155 + t315 + (-pkin(5) * t8 - t11 * t20 - t12 * t21 - t28 * t35) * m(7) + (m(7) * t283 + t269 + (-m(7) * t270 + t268) * qJD(6)) * pkin(11); -t425 + t424 - t28 * (mrSges(7,1) * t111 + mrSges(7,2) * t110) + (Ifges(7,1) * t110 - t361) * t389 + t46 * t388 + (Ifges(7,5) * t110 - Ifges(7,6) * t111) * t387 - t11 * t77 + t12 * t78 + (t11 * t110 + t111 * t12) * mrSges(7,3) + t13 + (-Ifges(7,2) * t111 + t106 + t47) * t391;];
tauc  = t3(:);

% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:16
% EndTime: 2019-03-09 07:12:46
% DurationCPUTime: 15.89s
% Computational Cost: add. (21533->661), mult. (55006->911), div. (0->0), fcn. (42743->10), ass. (0->286)
t272 = sin(pkin(11));
t273 = cos(pkin(11));
t277 = sin(qJ(3));
t281 = cos(qJ(3));
t251 = t272 * t281 + t273 * t277;
t242 = t251 * qJD(1);
t276 = sin(qJ(4));
t280 = cos(qJ(4));
t217 = qJD(3) * t280 - t242 * t276;
t218 = qJD(3) * t276 + t242 * t280;
t275 = sin(qJ(5));
t279 = cos(qJ(5));
t159 = t217 * t275 + t218 * t279;
t274 = sin(qJ(6));
t278 = cos(qJ(6));
t305 = t279 * t217 - t218 * t275;
t102 = t159 * t278 + t274 * t305;
t362 = pkin(7) + qJ(2);
t262 = t362 * t272;
t252 = qJD(1) * t262;
t263 = t362 * t273;
t253 = qJD(1) * t263;
t204 = -t252 * t281 - t277 * t253;
t197 = -qJD(3) * pkin(3) - t204;
t150 = -pkin(4) * t217 + t197;
t103 = -pkin(5) * t305 + t150;
t429 = pkin(10) * t305;
t436 = -t272 * t277 + t281 * t273;
t241 = t436 * qJD(1);
t311 = -pkin(2) * t273 - pkin(1);
t261 = qJD(1) * t311 + qJD(2);
t173 = -pkin(3) * t241 - pkin(8) * t242 + t261;
t205 = -t277 * t252 + t281 * t253;
t198 = qJD(3) * pkin(8) + t205;
t130 = t173 * t276 + t198 * t280;
t111 = pkin(9) * t217 + t130;
t106 = t279 * t111;
t129 = t280 * t173 - t198 * t276;
t110 = -pkin(9) * t218 + t129;
t234 = qJD(4) - t241;
t91 = pkin(4) * t234 + t110;
t56 = t275 * t91 + t106;
t41 = t56 + t429;
t345 = t274 * t41;
t229 = qJD(5) + t234;
t441 = pkin(10) * t159;
t104 = t275 * t111;
t55 = t279 * t91 - t104;
t40 = t55 - t441;
t38 = pkin(5) * t229 + t40;
t16 = t278 * t38 - t345;
t343 = t278 * t41;
t17 = t274 * t38 + t343;
t244 = t251 * qJD(3);
t231 = qJD(1) * t244;
t434 = -t159 * t274 + t278 * t305;
t243 = t436 * qJD(3);
t230 = qJD(1) * t243;
t167 = qJD(4) * t217 + t230 * t280;
t168 = -qJD(4) * t218 - t230 * t276;
t80 = qJD(5) * t305 + t167 * t279 + t168 * t275;
t81 = -qJD(5) * t159 - t167 * t275 + t168 * t279;
t32 = qJD(6) * t434 + t274 * t81 + t278 * t80;
t33 = -qJD(6) * t102 - t274 * t80 + t278 * t81;
t314 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t231;
t358 = Ifges(7,4) * t102;
t225 = qJD(6) + t229;
t377 = -t225 / 0.2e1;
t387 = t102 / 0.2e1;
t388 = -t102 / 0.2e1;
t392 = -t434 / 0.2e1;
t287 = t436 * qJD(2);
t160 = qJD(1) * t287 + qJD(3) * t204;
t181 = pkin(3) * t231 - pkin(8) * t230;
t75 = -qJD(4) * t130 - t160 * t276 + t280 * t181;
t53 = pkin(4) * t231 - pkin(9) * t167 + t75;
t319 = qJD(4) * t280;
t320 = qJD(4) * t276;
t74 = t280 * t160 + t173 * t319 + t276 * t181 - t198 * t320;
t60 = pkin(9) * t168 + t74;
t13 = -qJD(5) * t56 - t275 * t60 + t279 * t53;
t6 = pkin(5) * t231 - pkin(10) * t80 + t13;
t317 = qJD(5) * t279;
t318 = qJD(5) * t275;
t12 = -t111 * t318 + t275 * t53 + t279 * t60 + t91 * t317;
t9 = pkin(10) * t81 + t12;
t2 = qJD(6) * t16 + t274 * t6 + t278 * t9;
t3 = -qJD(6) * t17 - t274 * t9 + t278 * t6;
t428 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t51 = Ifges(7,2) * t434 + Ifges(7,6) * t225 + t358;
t96 = Ifges(7,4) * t434;
t52 = Ifges(7,1) * t102 + Ifges(7,5) * t225 + t96;
t458 = (Ifges(7,1) * t434 - t358) * t388 + (Ifges(7,5) * t434 - Ifges(7,6) * t102) * t377 + (t102 * t17 + t16 * t434) * mrSges(7,3) - t103 * (mrSges(7,1) * t102 + mrSges(7,2) * t434) + t51 * t387 + t314 + t428 + (-Ifges(7,2) * t102 + t52 + t96) * t392;
t199 = pkin(3) * t242 - pkin(8) * t241;
t141 = t280 * t199 - t204 * t276;
t390 = -pkin(9) - pkin(8);
t310 = qJD(4) * t390;
t367 = pkin(9) * t280;
t457 = -pkin(4) * t242 + t241 * t367 + t280 * t310 - t141;
t142 = t276 * t199 + t280 * t204;
t335 = t241 * t276;
t456 = -pkin(9) * t335 - t276 * t310 + t142;
t154 = Ifges(6,4) * t305;
t313 = Ifges(6,5) * t80 + Ifges(6,6) * t81 + Ifges(6,3) * t231;
t359 = Ifges(6,4) * t159;
t375 = -t229 / 0.2e1;
t384 = -t159 / 0.2e1;
t386 = -t305 / 0.2e1;
t417 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t89 = t159 * Ifges(6,1) + t229 * Ifges(6,5) + t154;
t455 = t313 + t417 + (Ifges(6,5) * t305 - Ifges(6,6) * t159) * t375 + (t159 * t56 + t305 * t55) * mrSges(6,3) + (-Ifges(6,2) * t159 + t154 + t89) * t386 - t150 * (mrSges(6,1) * t159 + mrSges(6,2) * t305) + (Ifges(6,1) * t305 - t359) * t384 + t458;
t255 = t275 * t280 + t276 * t279;
t179 = t255 * t241;
t405 = qJD(4) + qJD(5);
t210 = t405 * t255;
t410 = t179 - t210;
t291 = t275 * t276 - t279 * t280;
t180 = t291 * t241;
t209 = t405 * t291;
t409 = t180 - t209;
t264 = t390 * t276;
t265 = t390 * t280;
t216 = t275 * t264 - t279 * t265;
t416 = -qJD(5) * t216 + t456 * t275 + t279 * t457;
t415 = t264 * t317 + t265 * t318 + t275 * t457 - t456 * t279;
t446 = pkin(10) * t410 + t415;
t445 = -pkin(5) * t242 - pkin(10) * t409 + t416;
t353 = t218 * Ifges(5,4);
t144 = t217 * Ifges(5,2) + t234 * Ifges(5,6) + t353;
t214 = Ifges(5,4) * t217;
t145 = t218 * Ifges(5,1) + t234 * Ifges(5,5) + t214;
t294 = t129 * t280 + t130 * t276;
t297 = Ifges(5,5) * t280 - Ifges(5,6) * t276;
t360 = Ifges(5,4) * t280;
t299 = -Ifges(5,2) * t276 + t360;
t361 = Ifges(5,4) * t276;
t301 = Ifges(5,1) * t280 - t361;
t302 = mrSges(5,1) * t276 + mrSges(5,2) * t280;
t369 = t280 / 0.2e1;
t370 = -t276 / 0.2e1;
t378 = t218 / 0.2e1;
t444 = t145 * t369 + t144 * t370 + t197 * t302 + t217 * t299 / 0.2e1 + t301 * t378 + t234 * t297 / 0.2e1 - t294 * mrSges(5,3);
t268 = pkin(4) * t279 + pkin(5);
t315 = qJD(6) * t278;
t316 = qJD(6) * t274;
t328 = t275 * t278;
t61 = -t110 * t275 - t106;
t42 = t61 - t429;
t62 = t279 * t110 - t104;
t43 = t62 - t441;
t440 = t274 * t43 - t278 * t42 - t268 * t316 + (-t275 * t315 + (-t274 * t279 - t328) * qJD(5)) * pkin(4);
t329 = t274 * t275;
t439 = -t274 * t42 - t278 * t43 + t268 * t315 + (-t275 * t316 + (t278 * t279 - t329) * qJD(5)) * pkin(4);
t233 = Ifges(4,4) * t241;
t404 = t233 / 0.2e1 + t242 * Ifges(4,1) / 0.2e1;
t438 = t261 * mrSges(4,2) + Ifges(4,5) * qJD(3) + t404 + t444;
t437 = t243 * t276 + t251 * t319;
t400 = t32 / 0.2e1;
t399 = t33 / 0.2e1;
t394 = t80 / 0.2e1;
t393 = t81 / 0.2e1;
t88 = Ifges(6,2) * t305 + t229 * Ifges(6,6) + t359;
t431 = t88 / 0.2e1;
t373 = t231 / 0.2e1;
t430 = -t241 * Ifges(4,2) / 0.2e1;
t215 = t279 * t264 + t265 * t275;
t183 = -pkin(10) * t255 + t215;
t184 = -pkin(10) * t291 + t216;
t132 = t183 * t278 - t184 * t274;
t427 = qJD(6) * t132 + t274 * t445 + t278 * t446;
t133 = t183 * t274 + t184 * t278;
t426 = -qJD(6) * t133 - t274 * t446 + t278 * t445;
t192 = t291 * t251;
t107 = -mrSges(6,1) * t305 + mrSges(6,2) * t159;
t413 = m(6) * t150 + t107;
t206 = -t255 * t274 - t278 * t291;
t119 = qJD(6) * t206 - t209 * t278 - t210 * t274;
t125 = -t179 * t274 - t180 * t278;
t326 = t119 - t125;
t207 = t255 * t278 - t274 * t291;
t120 = -qJD(6) * t207 + t209 * t274 - t210 * t278;
t124 = -t179 * t278 + t180 * t274;
t325 = t120 - t124;
t163 = pkin(4) * t335 + t205;
t412 = pkin(4) * t320 - pkin(5) * t410 - t163;
t203 = -pkin(3) * t436 - pkin(8) * t251 + t311;
t213 = -t262 * t277 + t263 * t281;
t148 = t280 * t203 - t213 * t276;
t118 = -pkin(4) * t436 - t251 * t367 + t148;
t208 = t280 * t213;
t149 = t276 * t203 + t208;
t333 = t251 * t276;
t131 = -pkin(9) * t333 + t149;
t73 = t275 * t118 + t279 * t131;
t411 = Ifges(5,5) * t167 + Ifges(5,6) * t168;
t408 = -t281 * t262 - t263 * t277;
t407 = -t276 * t75 + t280 * t74;
t406 = t75 * mrSges(5,1) - t74 * mrSges(5,2);
t403 = (m(3) * qJ(2) + mrSges(3,3)) * (t272 ^ 2 + t273 ^ 2);
t402 = Ifges(7,4) * t400 + Ifges(7,2) * t399 + Ifges(7,6) * t373;
t401 = Ifges(7,1) * t400 + Ifges(7,4) * t399 + Ifges(7,5) * t373;
t398 = Ifges(6,4) * t394 + Ifges(6,2) * t393 + Ifges(6,6) * t373;
t397 = Ifges(6,1) * t394 + Ifges(6,4) * t393 + Ifges(6,5) * t373;
t391 = t434 / 0.2e1;
t385 = t305 / 0.2e1;
t383 = t159 / 0.2e1;
t382 = t167 / 0.2e1;
t381 = t168 / 0.2e1;
t380 = -t217 / 0.2e1;
t379 = -t218 / 0.2e1;
t376 = t225 / 0.2e1;
t374 = t229 / 0.2e1;
t372 = -t234 / 0.2e1;
t288 = t251 * qJD(2);
t161 = qJD(1) * t288 + qJD(3) * t205;
t338 = t161 * t408;
t324 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t217 + mrSges(5,2) * t218 + mrSges(4,3) * t242;
t312 = Ifges(5,3) * t231 + t411;
t269 = -pkin(4) * t280 - pkin(3);
t308 = t231 * mrSges(4,1) + t230 * mrSges(4,2);
t72 = t279 * t118 - t131 * t275;
t174 = qJD(3) * t408 + t287;
t200 = pkin(3) * t244 - pkin(8) * t243;
t306 = -t174 * t276 + t280 * t200;
t176 = pkin(4) * t333 - t408;
t303 = mrSges(5,1) * t280 - mrSges(5,2) * t276;
t300 = Ifges(5,1) * t276 + t360;
t298 = Ifges(5,2) * t280 + t361;
t296 = Ifges(5,5) * t276 + Ifges(5,6) * t280;
t57 = -pkin(5) * t436 + pkin(10) * t192 + t72;
t191 = t255 * t251;
t63 = -pkin(10) * t191 + t73;
t26 = -t274 * t63 + t278 * t57;
t27 = t274 * t57 + t278 * t63;
t295 = -t276 * t74 - t280 * t75;
t293 = t129 * t276 - t130 * t280;
t171 = -mrSges(5,2) * t234 + mrSges(5,3) * t217;
t172 = mrSges(5,1) * t234 - mrSges(5,3) * t218;
t292 = t171 * t280 - t172 * t276;
t135 = -t191 * t278 + t192 * t274;
t136 = -t191 * t274 - t192 * t278;
t66 = -t243 * t367 + pkin(4) * t244 + (-t208 + (pkin(9) * t251 - t203) * t276) * qJD(4) + t306;
t82 = t280 * t174 + t276 * t200 + t203 * t319 - t213 * t320;
t71 = -pkin(9) * t437 + t82;
t20 = t118 * t317 - t131 * t318 + t275 * t66 + t279 * t71;
t21 = -qJD(5) * t73 - t275 * t71 + t279 * t66;
t175 = qJD(3) * t213 + t288;
t138 = pkin(4) * t437 + t175;
t116 = -pkin(4) * t168 + t161;
t283 = t129 * mrSges(5,1) + t16 * mrSges(7,1) + t261 * mrSges(4,1) + t55 * mrSges(6,1) + t234 * Ifges(5,3) + t218 * Ifges(5,5) + t217 * Ifges(5,6) - Ifges(4,6) * qJD(3) - t242 * Ifges(4,4) + t430 + t225 * Ifges(7,3) + t102 * Ifges(7,5) + t434 * Ifges(7,6) + t229 * Ifges(6,3) + t159 * Ifges(6,5) + t305 * Ifges(6,6) - t130 * mrSges(5,2) - t17 * mrSges(7,2) - t56 * mrSges(6,2);
t237 = pkin(4) * t328 + t268 * t274;
t236 = -pkin(4) * t329 + t268 * t278;
t228 = pkin(5) * t291 + t269;
t219 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t241;
t147 = -mrSges(5,2) * t231 + mrSges(5,3) * t168;
t146 = mrSges(5,1) * t231 - mrSges(5,3) * t167;
t140 = mrSges(6,1) * t229 - mrSges(6,3) * t159;
t139 = -mrSges(6,2) * t229 + mrSges(6,3) * t305;
t137 = pkin(4) * t218 + pkin(5) * t159;
t134 = pkin(5) * t191 + t176;
t117 = -mrSges(5,1) * t168 + mrSges(5,2) * t167;
t95 = t167 * Ifges(5,1) + t168 * Ifges(5,4) + t231 * Ifges(5,5);
t94 = t167 * Ifges(5,4) + t168 * Ifges(5,2) + t231 * Ifges(5,6);
t93 = t192 * t405 - t255 * t243;
t92 = -t210 * t251 - t243 * t291;
t85 = mrSges(7,1) * t225 - mrSges(7,3) * t102;
t84 = -mrSges(7,2) * t225 + mrSges(7,3) * t434;
t83 = -qJD(4) * t149 + t306;
t77 = -mrSges(6,2) * t231 + mrSges(6,3) * t81;
t76 = mrSges(6,1) * t231 - mrSges(6,3) * t80;
t69 = -pkin(5) * t93 + t138;
t59 = -mrSges(7,1) * t434 + mrSges(7,2) * t102;
t47 = -pkin(5) * t81 + t116;
t39 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t37 = -qJD(6) * t136 - t274 * t92 + t278 * t93;
t36 = qJD(6) * t135 + t274 * t93 + t278 * t92;
t29 = -mrSges(7,2) * t231 + mrSges(7,3) * t33;
t28 = mrSges(7,1) * t231 - mrSges(7,3) * t32;
t19 = t278 * t40 - t345;
t18 = -t274 * t40 - t343;
t15 = pkin(10) * t93 + t20;
t14 = pkin(5) * t244 - pkin(10) * t92 + t21;
t10 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t5 = -qJD(6) * t27 + t14 * t278 - t15 * t274;
t4 = qJD(6) * t26 + t14 * t274 + t15 * t278;
t1 = [(-Ifges(6,5) * t192 + Ifges(7,5) * t136 - Ifges(6,6) * t191 + Ifges(7,6) * t135 + t251 * t297) * t373 - (mrSges(4,3) * t230 + t117) * t408 - (-mrSges(4,3) * t160 - Ifges(4,4) * t230 + Ifges(6,5) * t394 + Ifges(7,5) * t400 + Ifges(6,6) * t393 + Ifges(7,6) * t399 + (Ifges(6,3) + Ifges(7,3)) * t373 + t406 + t417 + t428) * t436 + (t430 + t283) * t244 + (Ifges(7,4) * t136 + Ifges(7,2) * t135) * t399 + (t95 * t369 + t94 * t370 + Ifges(4,1) * t230 - Ifges(4,4) * t231 + t301 * t382 + t299 * t381 + (mrSges(4,3) + t302) * t161 + t295 * mrSges(5,3) + (-t280 * t144 / 0.2e1 + t145 * t370 + t197 * t303 + t298 * t380 + t300 * t379 + t296 * t372 + t293 * mrSges(5,3)) * qJD(4)) * t251 + t324 * t175 + (t404 + t438) * t243 + (Ifges(7,1) * t136 + Ifges(7,4) * t135) * t400 + t93 * t431 + t311 * t308 - t192 * t397 - t191 * t398 + t136 * t401 + t135 * t402 + (Ifges(6,1) * t92 + Ifges(6,4) * t93) * t383 + (Ifges(6,4) * t92 + Ifges(6,2) * t93) * t385 + (Ifges(7,1) * t36 + Ifges(7,4) * t37) * t387 + (Ifges(7,4) * t36 + Ifges(7,2) * t37) * t391 + (Ifges(6,5) * t92 + Ifges(6,6) * t93) * t374 + (Ifges(7,5) * t36 + Ifges(7,6) * t37) * t376 + m(5) * (t129 * t83 + t130 * t82 + t148 * t75 + t149 * t74 + t175 * t197 - t338) + m(4) * (t160 * t213 + t174 * t205 - t175 * t204 - t338) - (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t231 * t436 - (t314 + t313 + t312 + t411) * t436 / 0.2e1 + 0.2e1 * t403 * qJD(2) * qJD(1) + (-Ifges(6,1) * t192 - Ifges(6,4) * t191) * t394 + (-Ifges(6,4) * t192 - Ifges(6,2) * t191) * t393 + (-t12 * t191 + t13 * t192 - t55 * t92 + t56 * t93) * mrSges(6,3) + t116 * (mrSges(6,1) * t191 - mrSges(6,2) * t192) + m(6) * (t116 * t176 + t12 * t73 + t13 * t72 + t138 * t150 + t20 * t56 + t21 * t55) + m(7) * (t103 * t69 + t134 * t47 + t16 * t5 + t17 * t4 + t2 * t27 + t26 * t3) + t72 * t76 + t73 * t77 + (-t204 * t243 - t205 * t244 - t213 * t231) * mrSges(4,3) + t69 * t59 + t37 * t51 / 0.2e1 + t36 * t52 / 0.2e1 + (t135 * t2 - t136 * t3 - t16 * t36 + t17 * t37) * mrSges(7,3) + t27 * t29 + t26 * t28 + t4 * t84 + t5 * t85 + t92 * t89 / 0.2e1 + t103 * (-mrSges(7,1) * t37 + mrSges(7,2) * t36) + t134 * t10 + t47 * (-mrSges(7,1) * t135 + mrSges(7,2) * t136) + t138 * t107 + t20 * t139 + t21 * t140 + t148 * t146 + t149 * t147 + t150 * (-mrSges(6,1) * t93 + mrSges(6,2) * t92) + t82 * t171 + t83 * t172 + t176 * t39 + t174 * t219; t308 + (-t219 - t292) * t241 + t292 * qJD(4) + t409 * t139 + t326 * t84 + t410 * t140 + t325 * t85 - m(4) * (-t204 * t242 + t205 * t241) + t206 * t28 + t207 * t29 + (-t59 - t107 - t324) * t242 - t291 * t76 + t255 * t77 + t276 * t147 + t280 * t146 - t403 * qJD(1) ^ 2 + (-t103 * t242 + t16 * t325 + t17 * t326 + t2 * t207 + t206 * t3) * m(7) + (t12 * t255 - t13 * t291 - t150 * t242 + t409 * t56 + t410 * t55) * m(6) + (-t197 * t242 - t234 * t293 - t295) * m(5); (-Ifges(6,5) * t180 - Ifges(6,6) * t179) * t375 + (-Ifges(6,1) * t180 - Ifges(6,4) * t179) * t384 + (-Ifges(6,4) * t180 - Ifges(6,2) * t179) * t386 + (-mrSges(7,1) * t325 + mrSges(7,2) * t326) * t103 + (-t16 * t326 + t17 * t325 + t2 * t206 - t207 * t3) * mrSges(7,3) + (t119 / 0.2e1 - t125 / 0.2e1) * t52 + (t120 / 0.2e1 - t124 / 0.2e1) * t51 + t415 * t139 + (t116 * t269 + t12 * t216 + t13 * t215 - t150 * t163 + t415 * t56 + t416 * t55) * m(6) + t416 * t140 + (-mrSges(4,1) - t303) * t161 - t324 * t205 + (t204 * mrSges(4,3) - t233 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t242 - t438) * t241 + t426 * t85 + t427 * t84 + (t103 * t412 + t132 * t3 + t133 * t2 + t16 * t426 + t17 * t427 + t228 * t47) * m(7) + (-pkin(3) * t161 - t129 * t141 - t130 * t142 - t197 * t205) * m(5) + t255 * t397 + (Ifges(7,4) * t207 + Ifges(7,2) * t206) * t399 + (Ifges(7,1) * t207 + Ifges(7,4) * t206) * t400 + t207 * t401 + t206 * t402 + t298 * t381 + t300 * t382 + (Ifges(7,1) * t119 + Ifges(7,4) * t120) * t387 + (Ifges(7,1) * t125 + Ifges(7,4) * t124) * t388 + (Ifges(7,4) * t119 + Ifges(7,2) * t120) * t391 + (Ifges(7,4) * t125 + Ifges(7,2) * t124) * t392 + t94 * t369 + (Ifges(7,5) * t119 + Ifges(7,6) * t120) * t376 + (Ifges(7,5) * t125 + Ifges(7,6) * t124) * t377 - t291 * t398 + (Ifges(6,5) * t255 + Ifges(7,5) * t207 - Ifges(6,6) * t291 + Ifges(7,6) * t206 + t296) * t373 + (Ifges(6,4) * t255 - Ifges(6,2) * t291) * t393 + (Ifges(6,1) * t255 - Ifges(6,4) * t291) * t394 + t116 * (mrSges(6,1) * t291 + mrSges(6,2) * t255) + (t180 / 0.2e1 - t209 / 0.2e1) * t89 + (-Ifges(6,1) * t209 - Ifges(6,4) * t210) * t383 + (-Ifges(6,4) * t209 - Ifges(6,2) * t210) * t385 + (-Ifges(6,5) * t209 - Ifges(6,6) * t210) * t374 + (t179 / 0.2e1 - t210 / 0.2e1) * t88 + ((-m(5) * t294 - t276 * t171 - t280 * t172) * qJD(4) - t146 * t276 + t147 * t280 + m(5) * t407) * pkin(8) + t407 * mrSges(5,3) + t412 * t59 + (-mrSges(6,1) * t410 + mrSges(6,2) * t409) * t150 + (-t12 * t291 - t13 * t255 - t409 * t55 + t410 * t56) * mrSges(6,3) + (t205 * mrSges(4,3) - t283) * t242 + (pkin(4) * t276 * t413 + t444) * qJD(4) - pkin(3) * t117 + t132 * t28 + t133 * t29 - t160 * mrSges(4,2) - t163 * t107 - t142 * t171 - t141 * t172 + t47 * (-mrSges(7,1) * t206 + mrSges(7,2) * t207) + t215 * t76 + t216 * t77 - t204 * t219 + t228 * t10 + Ifges(4,5) * t230 - Ifges(4,6) * t231 + t269 * t39 + t276 * t95 / 0.2e1; t406 + (-Ifges(5,2) * t218 + t145 + t214) * t380 + t455 + (Ifges(5,5) * t217 - Ifges(5,6) * t218) * t372 + t144 * t378 + (Ifges(5,1) * t217 - t353) * t379 - m(6) * (t55 * t61 + t56 * t62) + t159 * t431 + (t275 * t77 + t279 * t76 + (t139 * t279 - t140 * t275) * qJD(5) + m(6) * (t12 * t275 + t13 * t279 + t317 * t56 - t318 * t55) - t413 * t218) * pkin(4) + (t129 * t217 + t130 * t218) * mrSges(5,3) + t439 * t84 + t440 * t85 + (-t103 * t137 + t16 * t440 + t17 * t439 + t2 * t237 + t236 * t3) * m(7) + t312 - t137 * t59 - t62 * t139 - t61 * t140 - t129 * t171 + t130 * t172 - t197 * (mrSges(5,1) * t218 + mrSges(5,2) * t217) + t236 * t28 + t237 * t29; (-t159 * t59 + t274 * t29 + t278 * t28 + (-t274 * t85 + t278 * t84) * qJD(6) + (-t103 * t159 - t16 * t316 + t17 * t315 + t2 * t274 + t278 * t3) * m(7)) * pkin(5) - m(7) * (t16 * t18 + t17 * t19) + t88 * t383 - t19 * t84 - t18 * t85 - t55 * t139 + t56 * t140 + t455; -t16 * t84 + t17 * t85 + t458;];
tauc  = t1(:);

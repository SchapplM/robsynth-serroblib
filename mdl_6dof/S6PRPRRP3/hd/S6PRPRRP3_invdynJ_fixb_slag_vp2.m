% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:44
% EndTime: 2019-03-08 20:05:13
% DurationCPUTime: 16.57s
% Computational Cost: add. (6713->573), mult. (16113->751), div. (0->0), fcn. (12637->14), ass. (0->274)
t445 = Ifges(6,4) + Ifges(7,4);
t211 = sin(pkin(11));
t213 = cos(pkin(11));
t218 = sin(qJ(4));
t221 = cos(qJ(4));
t180 = t211 * t218 - t221 * t213;
t212 = sin(pkin(6));
t222 = cos(qJ(2));
t315 = t212 * t222;
t226 = t180 * t315;
t354 = pkin(8) + qJ(3);
t188 = t354 * t211;
t189 = t354 * t213;
t402 = -t221 * t188 - t189 * t218;
t413 = qJD(1) * t226 - t180 * qJD(3) + qJD(4) * t402;
t175 = t180 * qJD(4);
t181 = t211 * t221 + t213 * t218;
t176 = t181 * qJD(4);
t219 = sin(qJ(2));
t309 = qJD(1) * t212;
t284 = t219 * t309;
t447 = pkin(4) * t176 + pkin(9) * t175 - t284;
t220 = cos(qJ(5));
t205 = pkin(5) * t220 + pkin(4);
t446 = -m(6) * pkin(4) - m(7) * t205 - mrSges(5,1);
t215 = -qJ(6) - pkin(9);
t384 = -m(6) * pkin(9) + m(7) * t215 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t123 = -qJD(2) * t175 + qJDD(2) * t181;
t174 = t181 * qJD(2);
t217 = sin(qJ(5));
t143 = qJD(4) * t220 - t174 * t217;
t59 = qJD(5) * t143 + qJDD(4) * t217 + t123 * t220;
t375 = t59 / 0.2e1;
t144 = qJD(4) * t217 + t174 * t220;
t60 = -qJD(5) * t144 + qJDD(4) * t220 - t123 * t217;
t374 = t60 / 0.2e1;
t124 = -qJD(2) * t176 - qJDD(2) * t180;
t116 = qJDD(5) - t124;
t373 = t116 / 0.2e1;
t428 = Ifges(6,1) + Ifges(7,1);
t426 = Ifges(6,5) + Ifges(7,5);
t425 = Ifges(6,2) + Ifges(7,2);
t424 = Ifges(6,6) + Ifges(7,6);
t423 = Ifges(7,3) + Ifges(6,3);
t444 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t303 = qJD(5) * t217;
t173 = t180 * qJD(2);
t325 = t173 * t217;
t443 = t303 + t325;
t416 = mrSges(4,3) * (t211 ^ 2 + t213 ^ 2);
t442 = t426 * t373 + t445 * t374 + t428 * t375;
t441 = -t217 * t413 + t220 * t447;
t204 = pkin(3) * t213 + pkin(2);
t119 = pkin(4) * t180 - pkin(9) * t181 - t204;
t302 = qJD(5) * t220;
t440 = t119 * t302 + t217 * t447 + t220 * t413;
t439 = t445 * t143;
t438 = t445 * t144;
t437 = t445 * t220;
t436 = t445 * t217;
t435 = qJD(5) + t173;
t210 = pkin(11) + qJ(4);
t206 = sin(t210);
t207 = cos(t210);
t434 = t206 * t384 + t207 * t446;
t376 = m(7) * pkin(5);
t433 = t116 * t424 + t425 * t60 + t445 * t59;
t138 = -t188 * t218 + t189 * t221;
t132 = t220 * t138;
t237 = qJ(6) * t175 - qJD(6) * t181;
t432 = pkin(5) * t176 + t237 * t220 + (-t132 + (qJ(6) * t181 - t119) * t217) * qJD(5) + t441;
t282 = t181 * t302;
t431 = -qJ(6) * t282 + (-qJD(5) * t138 + t237) * t217 + t440;
t430 = -mrSges(7,1) - mrSges(6,1);
t429 = -mrSges(6,2) - mrSges(7,2);
t421 = t143 * t424 + t144 * t426 + t423 * t435;
t420 = t143 * t425 + t424 * t435 + t438;
t419 = t144 * t428 + t426 * t435 + t439;
t57 = t217 * t119 + t132;
t418 = -qJD(5) * t57 + t441;
t417 = -t138 * t303 + t440;
t415 = Ifges(5,4) * t173;
t414 = t173 * Ifges(5,2);
t348 = mrSges(5,3) * t174;
t412 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t143 - mrSges(6,2) * t144 - t348;
t267 = qJD(5) * t215;
t117 = pkin(4) * t174 + pkin(9) * t173;
t185 = qJD(2) * qJ(3) + t284;
t214 = cos(pkin(6));
t308 = qJD(1) * t214;
t200 = t213 * t308;
t338 = pkin(8) * qJD(2);
t133 = t200 + (-t185 - t338) * t211;
t147 = t213 * t185 + t211 * t308;
t134 = t213 * t338 + t147;
t65 = t133 * t221 - t218 * t134;
t31 = t217 * t117 + t220 * t65;
t411 = -qJ(6) * t325 + qJD(6) * t220 + t217 * t267 - t31;
t30 = t220 * t117 - t217 * t65;
t324 = t173 * t220;
t410 = -pkin(5) * t174 - qJ(6) * t324 - qJD(6) * t217 + t220 * t267 - t30;
t296 = qJDD(2) * t213;
t297 = qJDD(2) * t211;
t179 = -mrSges(4,1) * t296 + mrSges(4,2) * t297;
t44 = -t124 * mrSges(5,1) + t123 * mrSges(5,2);
t409 = t179 + t44;
t24 = -mrSges(6,1) * t60 + mrSges(6,2) * t59;
t408 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t123 + t24;
t407 = t376 + mrSges(7,1);
t66 = t133 * t218 + t134 * t221;
t406 = pkin(5) * t443 - t66;
t405 = Ifges(5,5) * qJD(4);
t404 = Ifges(5,6) * qJD(4);
t240 = t147 * t213 - (-t185 * t211 + t200) * t211;
t403 = t222 * t240;
t256 = -mrSges(4,1) * t213 + mrSges(4,2) * t211;
t401 = -mrSges(5,1) * t173 - mrSges(5,2) * t174 - t256 * qJD(2);
t251 = mrSges(7,1) * t217 + mrSges(7,2) * t220;
t253 = mrSges(6,1) * t217 + mrSges(6,2) * t220;
t61 = -qJD(4) * pkin(4) - t65;
t40 = -pkin(5) * t143 + qJD(6) + t61;
t400 = t40 * t251 + t61 * t253;
t399 = -t217 * t424 + t220 * t426;
t398 = -t217 * t425 + t437;
t397 = t220 * t428 - t436;
t396 = t116 * t423 + t424 * t60 + t426 * t59;
t395 = -t302 - t324;
t273 = qJD(2) * t309;
t195 = t222 * t273;
t300 = qJDD(1) * t212;
t163 = t219 * t300 + t195;
t145 = t163 + t444;
t299 = qJDD(1) * t214;
t198 = t213 * t299;
t112 = -t145 * t211 + t198;
t113 = t213 * t145 + t211 * t299;
t393 = -t112 * t211 + t113 * t213;
t304 = qJD(4) * t221;
t305 = qJD(4) * t218;
t93 = t198 + (-pkin(8) * qJDD(2) - t145) * t211;
t94 = pkin(8) * t296 + t113;
t16 = t133 * t304 - t134 * t305 + t218 * t93 + t221 * t94;
t14 = qJDD(4) * pkin(9) + t16;
t194 = t219 * t273;
t162 = t222 * t300 - t194;
t235 = qJDD(3) - t162;
t135 = -qJDD(2) * t204 + t235;
t39 = -pkin(4) * t124 - pkin(9) * t123 + t135;
t62 = qJD(4) * pkin(9) + t66;
t244 = -t222 * t309 + qJD(3);
t160 = -qJD(2) * t204 + t244;
t73 = pkin(4) * t173 - pkin(9) * t174 + t160;
t3 = t220 * t14 + t217 * t39 + t73 * t302 - t303 * t62;
t28 = t217 * t73 + t220 * t62;
t4 = -qJD(5) * t28 - t14 * t217 + t220 * t39;
t392 = -t217 * t4 + t220 * t3;
t389 = m(7) + m(6) + m(5);
t387 = -mrSges(6,1) - t407;
t275 = m(4) * qJ(3) + mrSges(4,3);
t386 = -mrSges(5,3) - t275 + mrSges(3,2);
t252 = -mrSges(7,1) * t220 + mrSges(7,2) * t217;
t254 = -mrSges(6,1) * t220 + mrSges(6,2) * t217;
t385 = -t252 - t254 - t446;
t74 = -mrSges(7,1) * t143 + mrSges(7,2) * t144;
t382 = -m(7) * t40 - t74;
t381 = -m(6) * t61 + t412;
t230 = m(4) * pkin(2) - t256;
t380 = t230 + mrSges(3,1) - t434;
t379 = t160 * mrSges(5,2) + t405 / 0.2e1;
t2 = qJ(6) * t60 + qJD(6) * t143 + t3;
t378 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t27 = -t217 * t62 + t220 * t73;
t20 = -qJ(6) * t144 + t27;
t12 = pkin(5) * t435 + t20;
t21 = qJ(6) * t143 + t28;
t377 = t21 * mrSges(7,2) + t28 * mrSges(6,2) + t404 / 0.2e1 - t12 * mrSges(7,1) - t160 * mrSges(5,1) - t27 * mrSges(6,1);
t372 = -t143 / 0.2e1;
t371 = t143 / 0.2e1;
t370 = -t144 / 0.2e1;
t369 = t144 / 0.2e1;
t368 = -t435 / 0.2e1;
t367 = t435 / 0.2e1;
t366 = t173 / 0.2e1;
t365 = -t174 / 0.2e1;
t364 = t174 / 0.2e1;
t358 = g(3) * t212;
t32 = mrSges(7,1) * t116 - mrSges(7,3) * t59;
t33 = mrSges(6,1) * t116 - mrSges(6,3) * t59;
t353 = t32 + t33;
t34 = -mrSges(7,2) * t116 + mrSges(7,3) * t60;
t35 = -mrSges(6,2) * t116 + mrSges(6,3) * t60;
t352 = t34 + t35;
t345 = mrSges(7,3) * t143;
t82 = -mrSges(7,2) * t435 + t345;
t347 = mrSges(6,3) * t143;
t83 = -mrSges(6,2) * t435 + t347;
t351 = t82 + t83;
t344 = mrSges(7,3) * t144;
t84 = mrSges(7,1) * t435 - t344;
t346 = mrSges(6,3) * t144;
t85 = mrSges(6,1) * t435 - t346;
t350 = t84 + t85;
t349 = mrSges(5,3) * t173;
t343 = Ifges(5,4) * t174;
t333 = cos(pkin(10));
t332 = sin(pkin(10));
t261 = t332 * t222;
t264 = t333 * t219;
t170 = t214 * t264 + t261;
t327 = t170 * t217;
t262 = t332 * t219;
t263 = t333 * t222;
t172 = -t214 * t262 + t263;
t326 = t172 * t217;
t323 = t175 * t217;
t322 = t175 * t220;
t321 = t181 * t217;
t320 = t181 * t220;
t318 = t207 * t217;
t317 = t207 * t220;
t316 = t212 * t219;
t314 = t217 * t222;
t313 = t220 * t222;
t307 = qJD(2) * t219;
t294 = -t74 + t412;
t290 = m(4) + t389;
t288 = t212 * t314;
t287 = t212 * t313;
t283 = t212 * t307;
t23 = -t60 * mrSges(7,1) + t59 * mrSges(7,2);
t269 = -t303 / 0.2e1;
t266 = t212 * t333;
t265 = t212 * t332;
t56 = t220 * t119 - t138 * t217;
t167 = -t211 * t316 + t213 * t214;
t168 = t211 * t214 + t213 * t316;
t239 = t221 * t167 - t168 * t218;
t98 = t167 * t218 + t168 * t221;
t17 = -t133 * t305 - t134 * t304 - t218 * t94 + t221 * t93;
t76 = -t217 * t98 - t287;
t236 = -t220 * t98 + t288;
t232 = t282 - t323;
t231 = t181 * t303 + t322;
t227 = t181 * t315;
t15 = -qJDD(4) * pkin(4) - t17;
t87 = qJD(3) * t181 + qJD(4) * t138;
t223 = qJD(2) ^ 2;
t191 = t215 * t220;
t190 = t215 * t217;
t183 = -qJD(2) * pkin(2) + t244;
t171 = t214 * t261 + t264;
t169 = -t214 * t263 + t262;
t156 = t206 * t214 + t207 * t316;
t155 = t206 * t316 - t214 * t207;
t153 = -qJD(4) * mrSges(5,2) - t349;
t152 = -qJDD(2) * pkin(2) + t235;
t128 = t172 * t207 + t206 * t265;
t127 = t172 * t206 - t207 * t265;
t126 = t170 * t207 - t206 * t266;
t125 = t170 * t206 + t207 * t266;
t103 = t174 * Ifges(5,1) + t405 - t415;
t102 = t343 + t404 - t414;
t101 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t124;
t88 = pkin(5) * t321 - t402;
t64 = qJD(2) * t227 + qJD(4) * t98;
t63 = -qJD(2) * t226 + qJD(4) * t239;
t43 = pkin(5) * t232 + t87;
t41 = -qJ(6) * t321 + t57;
t38 = pkin(5) * t180 - qJ(6) * t320 + t56;
t26 = qJD(5) * t236 - t217 * t63 + t220 * t283;
t25 = qJD(5) * t76 + t217 * t283 + t220 * t63;
t5 = -pkin(5) * t60 + qJDD(6) + t15;
t1 = pkin(5) * t116 - qJ(6) * t59 - qJD(6) * t144 + t4;
t6 = [t98 * t101 + t63 * t153 - t352 * t236 + t353 * t76 + t350 * t26 + t351 * t25 + (-t167 * t211 + t168 * t213) * qJDD(2) * mrSges(4,3) - (t23 + t408) * t239 - t294 * t64 + (-m(2) - m(3) - t290) * g(3) + m(5) * (t16 * t98 + t17 * t239 + t63 * t66 - t64 * t65) + m(4) * (t112 * t167 + t113 * t168) + m(7) * (t1 * t76 + t12 * t26 - t2 * t236 + t21 * t25 - t239 * t5 + t40 * t64) + m(6) * (-t15 * t239 - t236 * t3 + t25 * t28 + t26 * t27 + t4 * t76 + t61 * t64) + ((mrSges(3,1) * qJDD(2) - t409) * t222 + (-qJDD(2) * mrSges(3,2) - qJD(2) * t401) * t219 + m(5) * (-t135 * t222 + t160 * t307) + m(4) * (qJD(2) * t403 - t152 * t222 + t183 * t307) + m(3) * (t162 * t222 + t163 * t219) + (-t219 * mrSges(3,1) + (-mrSges(3,2) + t416) * t222) * t223) * t212 + (m(3) * t214 ^ 2 + m(2)) * qJDD(1); -t412 * t87 + t413 * t153 + (t231 * t27 - t232 * t28 - t3 * t321 - t320 * t4) * mrSges(6,3) + (-t282 / 0.2e1 + t323 / 0.2e1) * t420 + (-t135 * t204 + t138 * t16 - t160 * t284 + t17 * t402 + t413 * t66 - t65 * t87) * m(5) + (-t15 * t402 + t27 * t418 + t28 * t417 + t3 * t57 + t4 * t56 + t61 * t87) * m(6) - t408 * t402 + t401 * t284 + (-pkin(2) * t152 + t240 * qJD(3) + t393 * qJ(3) - (t183 * t219 + t403) * t309) * m(4) + (-t327 * t376 + t430 * (-t169 * t317 + t327) + t429 * (t169 * t318 + t170 * t220) - t389 * (-t169 * t204 + t170 * t354) + t386 * t170 + t380 * t169) * g(2) + (-t326 * t376 + t430 * (-t171 * t317 + t326) + t429 * (t171 * t318 + t172 * t220) - t389 * (-t171 * t204 + t172 * t354) + t386 * t172 + t380 * t171) * g(1) + (-t231 * t426 - t232 * t424) * t367 + t320 * t442 + (t65 * mrSges(5,3) - Ifges(5,1) * t364 + t415 / 0.2e1 - t103 / 0.2e1 - t379) * t175 + (-t195 + t444) * t416 + (m(5) * t65 + t381 + t382) * qJD(1) * t227 + t40 * (mrSges(7,1) * t232 - mrSges(7,2) * t231) + t61 * (mrSges(6,1) * t232 - mrSges(6,2) * t231) - (t219 * t275 + t222 * t230) * t358 + t393 * mrSges(4,3) + t417 * t83 + t418 * t85 - t419 * t322 / 0.2e1 + (-t222 * t358 + t162 + t194) * mrSges(3,1) + (-t1 * t320 + t12 * t231 - t2 * t321 - t21 * t232) * mrSges(7,3) + (t219 * t358 - t163 + t195) * mrSges(3,2) + (-t389 * t204 * t315 + (t430 * (t207 * t313 + t217 * t219) + t429 * (-t207 * t314 + t219 * t220) + t434 * t222 + (-t217 * t376 - t354 * t389 - mrSges(5,3)) * t219) * t212) * g(3) + t431 * t82 + t432 * t84 + (t1 * t38 + t12 * t432 + t2 * t41 + t21 * t431 + t40 * t43 + t5 * t88) * m(7) - t433 * t321 / 0.2e1 + t152 * t256 + (-t231 * t428 - t232 * t445) * t369 + (-t231 * t445 - t232 * t425) * t371 + (t135 * mrSges(5,2) - t17 * mrSges(5,3) + Ifges(5,1) * t123 + Ifges(5,4) * t124 + Ifges(5,5) * qJDD(4) + t15 * t253 + t251 * t5 + t269 * t419 + t373 * t399 + t374 * t398 + t375 * t397) * t181 + (-t66 * mrSges(5,3) - Ifges(5,4) * t364 + t414 / 0.2e1 - t102 / 0.2e1 + t424 * t371 + t426 * t369 + t423 * t367 - t377 + t421 / 0.2e1) * t176 + Ifges(3,3) * qJDD(2) + (Ifges(4,4) * t211 + Ifges(4,2) * t213) * t296 + (Ifges(4,1) * t211 + Ifges(4,4) * t213) * t297 - t204 * t44 - pkin(2) * t179 + t138 * t101 + t88 * t23 + t43 * t74 + t56 * t33 + t57 * t35 + t41 * t34 + t38 * t32 + (t396 / 0.2e1 + t135 * mrSges(5,1) + t1 * mrSges(7,1) - t16 * mrSges(5,3) - Ifges(5,4) * t123 - Ifges(5,2) * t124 - Ifges(5,6) * qJDD(4) + t423 * t373 + t424 * t374 + t426 * t375 + t378) * t180; t173 * t153 - t223 * t416 + t294 * t174 + (t351 * t435 + t353) * t220 + (-t350 * t435 + t352) * t217 + (t1 * t220 - t174 * t40 + t2 * t217 + t435 * (-t12 * t217 + t21 * t220)) * m(7) + (-t174 * t61 + t217 * t3 + t220 * t4 + t435 * (-t217 * t27 + t220 * t28)) * m(6) + (t173 * t66 + t174 * t65 + t135) * m(5) + (-qJD(2) * t240 + t152) * m(4) + (-g(1) * t171 - g(2) * t169 + g(3) * t315) * t290 + t409; t410 * t84 + (t1 * t190 + t12 * t410 - t191 * t2 - t205 * t5 + t21 * t411 + t40 * t406) * m(7) + t411 * t82 + (t302 / 0.2e1 + t324 / 0.2e1) * t419 + (t127 * t385 + t128 * t384) * g(1) + (t125 * t385 + t126 * t384) * g(2) + (t155 * t385 + t156 * t384) * g(3) + (-t349 - t153) * t65 + (-t415 + t103) * t366 + (t143 * t398 + t144 * t397 + t399 * t435) * qJD(5) / 0.2e1 - (Ifges(5,1) * t365 + t368 * t399 + t370 * t397 + t372 * t398 - t379 - t400) * t173 + t400 * qJD(5) + (-Ifges(5,2) * t366 + t368 * t423 + t370 * t426 + t372 * t424 + t377) * t174 + (t217 * t426 + t220 * t424) * t373 + t5 * t252 + t217 * t442 + (t348 + t381) * t66 + (-pkin(4) * t15 - t27 * t30 - t28 * t31) * m(6) + (m(6) * ((-t217 * t28 - t220 * t27) * qJD(5) + t392) - t85 * t302 - t83 * t303 - t217 * t33 + t220 * t35) * pkin(9) + (-t343 + t421) * t365 + (t217 * t428 + t437) * t375 + (t220 * t425 + t436) * t374 + t433 * t220 / 0.2e1 + t15 * t254 + (-t325 / 0.2e1 + t269) * t420 + t406 * t74 + (t27 * t395 - t28 * t443 + t392) * mrSges(6,3) + (-t1 * t217 + t12 * t395 + t2 * t220 - t21 * t443) * mrSges(7,3) + t102 * t364 + Ifges(5,3) * qJDD(4) - t205 * t23 + t190 * t32 - t191 * t34 + Ifges(5,6) * t124 + Ifges(5,5) * t123 - t30 * t85 - t31 * t83 - pkin(4) * t24 - t16 * mrSges(5,2) + t17 * mrSges(5,1); (t346 + t85) * t28 + t420 * t369 + (t429 * (-t126 * t220 - t169 * t217) + t387 * (-t126 * t217 + t169 * t220)) * g(2) + (t429 * (-t128 * t220 - t171 * t217) + t387 * (-t128 * t217 + t171 * t220)) * g(1) + (t429 * (-t156 * t220 + t288) + t387 * (-t156 * t217 - t287)) * g(3) + (t143 * t426 - t144 * t424) * t368 + (t143 * t428 - t438) * t370 + t378 + (t347 - t83) * t27 + (-t144 * t425 + t419 + t439) * t372 + (-m(7) * (-t12 + t20) + t344 + t84) * t21 + t407 * t1 + t12 * t345 - t61 * (mrSges(6,1) * t144 + mrSges(6,2) * t143) - t40 * (mrSges(7,1) * t144 + mrSges(7,2) * t143) - t20 * t82 + t396 + (t144 * t382 + t32) * pkin(5); -t143 * t82 + t144 * t84 + (-g(1) * t127 - g(2) * t125 - g(3) * t155 + t12 * t144 - t21 * t143 + t5) * m(7) + t23;];
tau  = t6;

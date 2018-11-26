% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:19:37
% EndTime: 2018-11-23 17:19:55
% DurationCPUTime: 18.09s
% Computational Cost: add. (10355->755), mult. (26685->1009), div. (0->0), fcn. (19174->8), ass. (0->326)
t445 = Ifges(7,4) + Ifges(6,4);
t446 = Ifges(7,1) + Ifges(6,1);
t444 = Ifges(7,5) + Ifges(6,5);
t443 = Ifges(7,2) + Ifges(6,2);
t442 = Ifges(7,6) + Ifges(6,6);
t259 = cos(pkin(6));
t345 = qJD(1) * t259;
t250 = qJD(2) + t345;
t264 = cos(qJ(4));
t261 = sin(qJ(4));
t265 = cos(qJ(2));
t258 = sin(pkin(6));
t346 = qJD(1) * t258;
t320 = t265 * t346;
t308 = t261 * t320;
t195 = t250 * t264 - t308;
t262 = sin(qJ(2));
t321 = t262 * t346;
t235 = qJD(4) + t321;
t260 = sin(qJ(5));
t263 = cos(qJ(5));
t131 = -t195 * t260 + t235 * t263;
t194 = -t250 * t261 - t264 * t320;
t191 = qJD(5) - t194;
t132 = t195 * t263 + t235 * t260;
t464 = t445 * t132;
t437 = t443 * t131 + t442 * t191 + t464;
t466 = -t437 / 0.2e1;
t465 = t445 * t131;
t330 = pkin(1) * t345;
t247 = t265 * t330;
t266 = -pkin(2) - pkin(9);
t411 = pkin(3) + pkin(8);
t124 = t250 * t266 + t321 * t411 + qJD(3) - t247;
t313 = -qJ(3) * t262 - pkin(1);
t172 = (t265 * t266 + t313) * t258;
t152 = qJD(1) * t172;
t80 = t124 * t261 + t152 * t264;
t63 = pkin(10) * t235 + t80;
t208 = pkin(8) * t320 + t262 * t330;
t180 = pkin(3) * t320 + t208;
t240 = t250 * qJ(3);
t145 = t240 + t180;
t83 = -pkin(4) * t194 - pkin(10) * t195 + t145;
t24 = -t260 * t63 + t263 * t83;
t17 = -qJ(6) * t132 + t24;
t10 = pkin(5) * t191 + t17;
t25 = t260 * t83 + t263 * t63;
t18 = qJ(6) * t131 + t25;
t283 = t24 * t263 + t25 * t260;
t298 = mrSges(7,1) * t260 + mrSges(7,2) * t263;
t300 = mrSges(6,1) * t260 + mrSges(6,2) * t263;
t385 = t263 / 0.2e1;
t388 = -t260 / 0.2e1;
t79 = t124 * t264 - t261 * t152;
t62 = -pkin(4) * t235 - t79;
t42 = -pkin(5) * t131 + qJD(6) + t62;
t436 = t446 * t132 + t444 * t191 + t465;
t463 = -t42 * t298 - t62 * t300 - t385 * t436 - t388 * t437 + t283 * mrSges(6,3) + (t10 * t263 + t18 * t260) * mrSges(7,3);
t462 = -t194 * Ifges(5,2) / 0.2e1;
t350 = t262 * t263;
t188 = (t260 * t265 + t261 * t350) * t346;
t226 = pkin(4) * t261 - pkin(10) * t264 + qJ(3);
t351 = t261 * t266;
t193 = t260 * t226 + t263 * t351;
t304 = pkin(4) * t264 + pkin(10) * t261;
t222 = qJD(4) * t304 + qJD(3);
t274 = -qJD(5) * t193 + t263 * t222;
t309 = t264 * t321;
t312 = -t260 * t266 + pkin(5);
t334 = qJD(6) * t263;
t336 = qJD(5) * t260;
t339 = qJD(4) * t263;
t119 = t247 + (-t304 - t411) * t321;
t242 = pkin(2) * t321;
t285 = pkin(9) * t262 - qJ(3) * t265;
t178 = t285 * t346 + t242;
t107 = t264 * t178 + t261 * t180;
t95 = pkin(10) * t320 + t107;
t49 = t263 * t119 - t260 * t95;
t461 = pkin(5) * t309 - t49 + (qJD(4) * t312 - t334) * t264 + t274 + (t261 * t339 + t264 * t336 + t188) * qJ(6);
t352 = t260 * t262;
t187 = (-t261 * t352 + t263 * t265) * t346;
t335 = qJD(5) * t263;
t315 = t264 * t335;
t337 = qJD(4) * t266;
t316 = t264 * t337;
t323 = t260 * t222 + t226 * t335 + t263 * t316;
t50 = t260 * t119 + t263 * t95;
t460 = -qJ(6) * t187 - t50 - qJ(6) * t315 + (-qJD(6) * t264 + (qJ(6) * qJD(4) - qJD(5) * t266) * t261) * t260 + t323;
t340 = qJD(4) * t261;
t106 = -t261 * t178 + t180 * t264;
t94 = -pkin(4) * t320 - t106;
t459 = t261 * t337 - t94 + (-t260 * t340 + t187 + t315) * pkin(5);
t392 = t191 / 0.2e1;
t401 = t132 / 0.2e1;
t403 = t131 / 0.2e1;
t456 = t445 * t260;
t426 = t263 * t446 - t456;
t457 = t445 * t263;
t428 = -t260 * t443 + t457;
t430 = -t260 * t442 + t263 * t444;
t458 = t392 * t430 + t401 * t426 + t403 * t428 - t463;
t207 = pkin(8) * t321 - t247;
t455 = -qJD(3) - t207;
t454 = -t235 * Ifges(5,6) / 0.2e1 - Ifges(5,4) * t195 / 0.2e1;
t453 = t462 + t454;
t449 = -t250 / 0.2e1;
t314 = qJD(2) * t346;
t307 = t262 * t314;
t338 = qJD(4) * t264;
t141 = -qJD(4) * t308 + t250 * t338 - t264 * t307;
t344 = qJD(2) * t262;
t317 = t261 * t344;
t140 = -t250 * t340 + (-t265 * t338 + t317) * t346;
t306 = t265 * t314;
t73 = qJD(5) * t131 + t140 * t263 + t260 * t306;
t74 = -qJD(5) * t132 - t140 * t260 + t263 * t306;
t438 = t141 * t444 + t445 * t74 + t446 * t73;
t452 = t438 / 0.2e1;
t451 = mrSges(4,1) + mrSges(3,3);
t450 = mrSges(4,2) - mrSges(3,1);
t439 = t141 * t442 + t443 * t74 + t445 * t73;
t448 = t321 / 0.2e1;
t447 = -t346 / 0.2e1;
t441 = Ifges(6,3) + Ifges(7,3);
t11 = Ifges(7,5) * t73 + Ifges(7,6) * t74 + Ifges(7,3) * t141;
t12 = Ifges(6,5) * t73 + Ifges(6,6) * t74 + Ifges(6,3) * t141;
t440 = t12 + t11;
t380 = -qJ(6) - pkin(10);
t311 = qJD(5) * t380;
t123 = pkin(4) * t195 - pkin(10) * t194;
t40 = t263 * t123 - t260 * t79;
t435 = -pkin(5) * t195 - qJD(6) * t260 - t40 + (qJ(6) * t194 + t311) * t263;
t355 = t194 * t260;
t41 = t260 * t123 + t263 * t79;
t434 = qJ(6) * t355 + t260 * t311 + t334 - t41;
t254 = t259 * t262 * pkin(1);
t353 = t258 * t265;
t218 = pkin(8) * t353 + t254;
t198 = -t259 * qJ(3) - t218;
t171 = pkin(3) * t353 - t198;
t215 = t259 * t261 + t264 * t353;
t325 = t261 * t353;
t216 = t259 * t264 - t325;
t105 = pkin(4) * t215 - pkin(10) * t216 + t171;
t354 = t258 * t262;
t251 = pkin(8) * t354;
t384 = pkin(1) * t265;
t322 = -pkin(2) - t384;
t155 = pkin(3) * t354 + t251 + (-pkin(9) + t322) * t259;
t98 = t261 * t155 + t264 * t172;
t92 = pkin(10) * t354 + t98;
t39 = t260 * t105 + t263 * t92;
t189 = Ifges(5,4) * t194;
t364 = t235 * Ifges(5,5);
t366 = t195 * Ifges(5,1);
t110 = t189 + t364 + t366;
t433 = t79 * mrSges(5,3) - t110 / 0.2e1 - t145 * mrSges(5,2) - t364 / 0.2e1 - t189 / 0.2e1;
t431 = t260 * t444 + t263 * t442;
t429 = t263 * t443 + t456;
t427 = t260 * t446 + t457;
t233 = pkin(2) * t307;
t342 = qJD(3) * t262;
t273 = (qJD(2) * t285 - t342) * t258;
t129 = qJD(1) * t273 + t233;
t181 = (t353 * t411 + t254) * qJD(2);
t156 = qJD(1) * t181;
t31 = t124 * t338 + t264 * t129 - t152 * t340 + t261 * t156;
t27 = pkin(10) * t306 + t31;
t234 = qJD(2) * t247;
t239 = t250 * qJD(3);
t310 = t411 * t354;
t282 = qJD(2) * t310;
t130 = -qJD(1) * t282 + t234 + t239;
t53 = pkin(4) * t141 - pkin(10) * t140 + t130;
t5 = t260 * t53 + t263 * t27 + t83 * t335 - t336 * t63;
t6 = -qJD(5) * t25 - t260 * t27 + t263 * t53;
t302 = -t260 * t6 + t263 * t5;
t32 = -t124 * t340 - t261 * t129 - t152 * t338 + t156 * t264;
t425 = t32 * mrSges(5,1) - t31 * mrSges(5,2) + Ifges(5,5) * t140 - Ifges(5,6) * t141;
t424 = -Ifges(3,4) * t320 / 0.2e1 + Ifges(3,5) * t449;
t326 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t327 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t329 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t55 = t132 * Ifges(7,5) + t131 * Ifges(7,6) + t191 * Ifges(7,3);
t56 = t132 * Ifges(6,5) + t131 * Ifges(6,6) + t191 * Ifges(6,3);
t422 = t327 * t131 - t329 * t132 + t326 * t191 - t18 * mrSges(7,2) - t25 * mrSges(6,2) - t80 * mrSges(5,3) + t453 + t55 / 0.2e1 + t56 / 0.2e1 + t10 * mrSges(7,1) + t145 * mrSges(5,1) + t24 * mrSges(6,1) + t454;
t421 = m(5) / 0.2e1;
t420 = m(6) / 0.2e1;
t419 = m(7) / 0.2e1;
t418 = Ifges(3,5) / 0.2e1;
t417 = Ifges(4,5) / 0.2e1;
t416 = Ifges(5,2) / 0.2e1;
t415 = t73 / 0.2e1;
t414 = t74 / 0.2e1;
t410 = m(5) * t79;
t409 = m(5) * t80;
t408 = m(6) * t62;
t407 = m(7) * t42;
t406 = pkin(1) * mrSges(3,1);
t405 = pkin(1) * mrSges(3,2);
t404 = -t131 / 0.2e1;
t402 = -t132 / 0.2e1;
t400 = t141 / 0.2e1;
t393 = -t191 / 0.2e1;
t391 = -t215 / 0.2e1;
t389 = t216 / 0.2e1;
t383 = pkin(5) * t260;
t374 = Ifges(4,6) * t265;
t371 = t140 * Ifges(5,1);
t370 = t140 * Ifges(5,4);
t369 = t141 * Ifges(5,4);
t360 = t250 * Ifges(4,5);
t120 = mrSges(5,1) * t306 - mrSges(5,3) * t140;
t22 = -mrSges(6,1) * t74 + mrSges(6,2) * t73;
t359 = t120 - t22;
t147 = mrSges(5,1) * t235 - mrSges(5,3) * t195;
t85 = -mrSges(6,1) * t131 + mrSges(6,2) * t132;
t358 = t85 - t147;
t100 = -mrSges(6,2) * t191 + mrSges(6,3) * t131;
t99 = -mrSges(7,2) * t191 + mrSges(7,3) * t131;
t357 = -t99 - t100;
t356 = qJ(6) * t264;
t101 = mrSges(7,1) * t191 - mrSges(7,3) * t132;
t102 = mrSges(6,1) * t191 - mrSges(6,3) * t132;
t349 = -t101 - t102;
t122 = -mrSges(5,1) * t194 + mrSges(5,2) * t195;
t203 = -mrSges(4,1) * t320 - mrSges(4,3) * t250;
t348 = t122 - t203;
t347 = t250 * t450 + t321 * t451;
t343 = qJD(2) * t265;
t341 = qJD(4) * t260;
t333 = t259 * t384;
t332 = -qJD(2) + t250 / 0.2e1;
t328 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t324 = t260 * t351;
t319 = t258 * t344;
t318 = t258 * t343;
t21 = -t74 * mrSges(7,1) + t73 * mrSges(7,2);
t38 = t263 * t105 - t260 * t92;
t97 = t155 * t264 - t261 * t172;
t1 = pkin(5) * t141 - qJ(6) * t73 - qJD(6) * t132 + t6;
t2 = qJ(6) * t74 + qJD(6) * t131 + t5;
t303 = -t1 * t260 + t2 * t263;
t301 = mrSges(6,1) * t263 - mrSges(6,2) * t260;
t299 = mrSges(7,1) * t263 - mrSges(7,2) * t260;
t248 = qJD(2) * t333;
t209 = -pkin(8) * t319 + t248;
t244 = pkin(2) * t319;
t149 = t244 + t273;
t44 = -t261 * t149 - t155 * t340 - t172 * t338 + t181 * t264;
t168 = -t216 * t260 + t258 * t350;
t169 = t216 * t263 + t258 * t352;
t43 = t264 * t149 + t155 * t338 - t172 * t340 + t261 * t181;
t36 = pkin(10) * t318 + t43;
t256 = t259 * qJD(3);
t154 = t248 + t256 - t282;
t166 = -qJD(4) * t215 + t258 * t317;
t167 = -qJD(4) * t325 + t259 * t338 - t264 * t319;
t72 = pkin(4) * t167 - pkin(10) * t166 + t154;
t7 = t105 * t335 + t260 * t72 + t263 * t36 - t336 * t92;
t196 = -pkin(8) * t307 + t234;
t91 = -pkin(4) * t354 - t97;
t199 = (-pkin(2) * t265 + t313) * t258;
t279 = (-qJ(3) * t343 - t342) * t258;
t184 = qJD(1) * t199;
t278 = Ifges(3,6) * t449 + (Ifges(3,4) * t262 + Ifges(3,2) * t265) * t447 + t360 / 0.2e1 + (-Ifges(4,6) * t262 - Ifges(4,3) * t265) * t346 / 0.2e1 - t184 * mrSges(4,2) - t208 * mrSges(3,3);
t210 = t218 * qJD(2);
t277 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t162 = -t196 - t239;
t197 = qJD(1) * t210;
t276 = -t196 * mrSges(3,2) - t162 * mrSges(4,3) + t197 * t450;
t37 = -pkin(4) * t318 - t44;
t8 = -qJD(5) * t39 - t260 * t36 + t263 * t72;
t28 = -pkin(4) * t306 - t32;
t271 = -t366 / 0.2e1 + t433;
t270 = t207 * mrSges(3,3) + t79 * mrSges(5,1) + t235 * Ifges(5,3) + t195 * Ifges(5,5) + t194 * Ifges(5,6) + Ifges(3,1) * t448 + Ifges(4,4) * t449 + (-Ifges(4,2) * t262 - t374) * t447 - t184 * mrSges(4,3) - t80 * mrSges(5,2) - t424;
t268 = t462 + t422;
t255 = -pkin(5) * t263 - pkin(4);
t238 = t380 * t263;
t237 = t380 * t260;
t232 = Ifges(3,5) * t306;
t231 = Ifges(4,5) * t307;
t230 = Ifges(5,3) * t306;
t223 = (-t266 + t383) * t264;
t221 = t263 * t226;
t217 = -t251 + t333;
t206 = -qJ(3) * t320 + t242;
t205 = (mrSges(4,2) * t265 - mrSges(4,3) * t262) * t346;
t202 = -mrSges(3,2) * t250 + mrSges(3,3) * t320;
t200 = t259 * t322 + t251;
t192 = t221 - t324;
t190 = -t209 - t256;
t186 = t250 * t260 - t263 * t309;
t185 = t250 * t263 + t260 * t309;
t182 = t244 + t279;
t179 = -qJD(1) * t310 + t247;
t177 = -t240 - t208;
t170 = -pkin(2) * t250 - t455;
t163 = -t260 * t356 + t193;
t159 = qJD(1) * t279 + t233;
t150 = t261 * t312 - t263 * t356 + t221;
t146 = -mrSges(5,2) * t235 + mrSges(5,3) * t194;
t121 = -mrSges(5,2) * t306 - mrSges(5,3) * t141;
t114 = -t260 * t316 + t274;
t113 = -qJD(5) * t324 + t323;
t90 = -qJD(5) * t169 - t166 * t260 + t263 * t318;
t89 = qJD(5) * t168 + t166 * t263 + t260 * t318;
t87 = mrSges(5,1) * t141 + mrSges(5,2) * t140;
t84 = -mrSges(7,1) * t131 + mrSges(7,2) * t132;
t77 = Ifges(5,5) * t306 - t369 + t371;
t76 = -t141 * Ifges(5,2) + Ifges(5,6) * t306 + t370;
t61 = -pkin(5) * t168 + t91;
t54 = pkin(5) * t355 + t80;
t48 = -mrSges(6,2) * t141 + mrSges(6,3) * t74;
t47 = -mrSges(7,2) * t141 + mrSges(7,3) * t74;
t46 = mrSges(6,1) * t141 - mrSges(6,3) * t73;
t45 = mrSges(7,1) * t141 - mrSges(7,3) * t73;
t29 = qJ(6) * t168 + t39;
t20 = pkin(5) * t215 - qJ(6) * t169 + t38;
t19 = -pkin(5) * t90 + t37;
t9 = -pkin(5) * t74 + t28;
t4 = qJ(6) * t90 + qJD(6) * t168 + t7;
t3 = pkin(5) * t167 - qJ(6) * t89 - qJD(6) * t169 + t8;
t13 = [(t56 + t55) * t167 / 0.2e1 + ((Ifges(5,5) * t389 + Ifges(5,6) * t391 - t199 * mrSges(4,3) + t200 * mrSges(4,1) - t217 * mrSges(3,3) + (-Ifges(4,4) + t418) * t259 + (t265 * t328 - 0.2e1 * t405) * t258) * t265 + (t198 * mrSges(4,1) - t199 * mrSges(4,2) - t218 * mrSges(3,3) + (t417 - Ifges(3,6)) * t259 + (-t262 * t328 - 0.2e1 * t406) * t258 + (Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t353) * t262) * t314 + (t167 * t444 + t445 * t90 + t446 * t89) * t401 + (t168 * t445 + t169 * t446 + t215 * t444) * t415 + t235 * (Ifges(5,5) * t166 - Ifges(5,6) * t167) / 0.2e1 + t194 * (Ifges(5,4) * t166 - Ifges(5,2) * t167) / 0.2e1 - t141 * (Ifges(5,4) * t216 - Ifges(5,2) * t215) / 0.2e1 + t436 * t89 / 0.2e1 + t437 * t90 / 0.2e1 + t140 * (Ifges(5,1) * t216 - Ifges(5,4) * t215) / 0.2e1 + t195 * (Ifges(5,1) * t166 - Ifges(5,4) * t167) / 0.2e1 + t439 * t168 / 0.2e1 + t440 * t215 / 0.2e1 + (t167 * t441 + t442 * t90 + t444 * t89) * t392 + (t168 * t442 + t169 * t444 + t215 * t441) * t400 + (t167 * t442 + t443 * t90 + t445 * t89) * t403 + (t168 * t443 + t169 * t445 + t215 * t442) * t414 + t347 * t210 + t167 * t453 + t169 * t452 + t77 * t389 + ((-mrSges(4,1) * t162 + mrSges(4,2) * t159 + mrSges(3,3) * t196) * t265 + (t230 / 0.2e1 - t159 * mrSges(4,3) + t451 * t197 + t425) * t262 + ((t177 * mrSges(4,1) + (t417 - Ifges(3,6) / 0.2e1) * t250 + t278) * t262 + (t170 * mrSges(4,1) + (-Ifges(4,4) / 0.2e1 + t418) * t250 + t270) * t265) * qJD(2)) * t258 + (t232 / 0.2e1 + t231 / 0.2e1 + t276) * t259 + (-t166 * t79 - t167 * t80 - t215 * t31 - t216 * t32) * mrSges(5,3) + t76 * t391 + t20 * t45 + t38 * t46 + t29 * t47 + t39 * t48 + t61 * t21 + t19 * t84 + t37 * t85 + t42 * (-mrSges(7,1) * t90 + mrSges(7,2) * t89) + t62 * (-mrSges(6,1) * t90 + mrSges(6,2) * t89) + t91 * t22 + t4 * t99 + t7 * t100 + t3 * t101 + t8 * t102 + t97 * t120 + t98 * t121 + t43 * t146 + t44 * t147 + t154 * t122 + t166 * t110 / 0.2e1 + t10 * (mrSges(7,1) * t167 - mrSges(7,3) * t89) + t24 * (mrSges(6,1) * t167 - mrSges(6,3) * t89) + t18 * (-mrSges(7,2) * t167 + mrSges(7,3) * t90) + t25 * (-mrSges(6,2) * t167 + mrSges(6,3) * t90) + t145 * (mrSges(5,1) * t167 + mrSges(5,2) * t166) + t9 * (-mrSges(7,1) * t168 + mrSges(7,2) * t169) + t28 * (-mrSges(6,1) * t168 + mrSges(6,2) * t169) + t171 * t87 + m(7) * (t1 * t20 + t10 * t3 + t18 * t4 + t19 * t42 + t2 * t29 + t61 * t9) + m(6) * (t24 * t8 + t25 * t7 + t28 * t91 + t37 * t62 + t38 * t6 + t39 * t5) + m(5) * (t130 * t171 + t145 * t154 + t31 * t98 + t32 * t97 + t43 * t80 + t44 * t79) + m(4) * (t159 * t199 + t162 * t198 + t170 * t210 + t177 * t190 + t182 * t184 + t197 * t200) + m(3) * (t196 * t218 - t197 * t217 + t207 * t210 + t208 * t209) + t190 * t203 + t182 * t205 + t209 * t202 + t2 * (-mrSges(7,2) * t215 + mrSges(7,3) * t168) + t5 * (-mrSges(6,2) * t215 + mrSges(6,3) * t168) + t1 * (mrSges(7,1) * t215 - mrSges(7,3) * t169) + t6 * (mrSges(6,1) * t215 - mrSges(6,3) * t169) + t130 * (mrSges(5,1) * t215 + mrSges(5,2) * t216); (t1 * t150 + t10 * t461 + t163 * t2 + t18 * t460 + t223 * t9 + t42 * t459) * m(7) + (((t146 + t409) * t266 + t268) * t264 + (t271 + (t358 + t408 - t410) * t266 + t428 * t404 + t426 * t402 + t430 * t393 + t463) * t261) * qJD(4) + ((-t270 + (t405 - t374 / 0.2e1) * t346 + qJD(2) * (Ifges(5,5) * t264 - Ifges(5,6) * t261) / 0.2e1 + (-pkin(2) * qJD(2) - t170) * mrSges(4,1) + t332 * Ifges(4,4) + t424) * t265 + (-t360 / 0.2e1 + ((t406 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t262) * t258 + (-Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t353) * qJD(1) + t332 * Ifges(3,6) + (-qJ(3) * qJD(2) - t177) * mrSges(4,1) + t271 * t261 + t268 * t264 - t278) * t262) * t346 + (-t32 * mrSges(5,3) + t371 / 0.2e1 - t369 / 0.2e1 + t130 * mrSges(5,2) + t9 * t298 + t28 * t300 + t77 / 0.2e1 + (-t1 * t263 - t2 * t260) * mrSges(7,3) + (-t260 * t5 - t263 * t6) * mrSges(6,3) + (m(5) * t32 - m(6) * t28 + t359) * t266 + (t42 * t299 + t62 * t301 + (t10 * t260 - t18 * t263) * mrSges(7,3) + (t24 * t260 - t25 * t263) * mrSges(6,3) + t429 * t404 + t427 * t402 + t431 * t393 + t263 * t466) * qJD(5) + t426 * t415 + t428 * t414 + t430 * t400 + (qJD(5) * t436 + t439) * t388 + t438 * t385) * t264 + (-t203 + t202) * t207 + t231 + t232 + (-t50 + t113) * t100 + (t130 * mrSges(5,1) - t370 / 0.2e1 + t11 / 0.2e1 + t12 / 0.2e1 - t76 / 0.2e1 - t31 * mrSges(5,3) + t327 * t74 - t329 * t73 + (m(5) * t31 + t121) * t266 + (t416 + t326) * t141 + t277) * t261 + t276 + (-pkin(2) * t197 - qJ(3) * t162 - t170 * t208 + t177 * t455 - t184 * t206) * m(4) - t436 * t188 / 0.2e1 + t460 * t99 + t461 * t101 - m(6) * (t24 * t49 + t25 * t50 + t62 * t94) - m(5) * (t106 * t79 + t107 * t80 + t145 * t179) + (-t49 + t114) * t102 + t459 * t84 + (t10 * t188 - t18 * t187) * mrSges(7,3) + t348 * qJD(3) - t347 * t208 + m(6) * (t25 * t113 + t24 * t114 + t6 * t192 + t5 * t193) + (Ifges(7,5) * t188 + Ifges(7,6) * t187) * t393 + (Ifges(6,5) * t188 + Ifges(6,6) * t187) * t393 + (Ifges(7,1) * t188 + Ifges(7,4) * t187) * t402 + (Ifges(6,1) * t188 + Ifges(6,4) * t187) * t402 + (Ifges(7,4) * t188 + Ifges(7,2) * t187) * t404 + (Ifges(6,4) * t188 + Ifges(6,2) * t187) * t404 + (-t187 * t25 + t188 * t24) * mrSges(6,3) + m(5) * (t130 * qJ(3) + t145 * qJD(3)) + qJ(3) * t87 - t94 * t85 - t107 * t146 - t106 * t147 + t150 * t45 + t163 * t47 - t179 * t122 + t187 * t466 - t42 * (-mrSges(7,1) * t187 + mrSges(7,2) * t188) - t62 * (-mrSges(6,1) * t187 + mrSges(6,2) * t188) + t192 * t46 + t193 * t48 - t206 * t205 + t223 * t21; t357 * t186 + t349 * t185 + (mrSges(4,1) * t343 + t205 * t262) * t346 + (t146 * t321 - t21 + (t260 * t349 - t263 * t357 + t146) * qJD(4) + t359) * t264 + (t121 + (t47 + t48) * t263 + (-t45 - t46) * t260 + (t260 * t357 + t263 * t349) * qJD(5) + t235 * (t84 + t358)) * t261 - m(7) * (t10 * t185 + t18 * t186) - m(6) * (t185 * t24 + t186 * t25) + 0.2e1 * ((qJD(4) * t80 + t32) * t421 + t409 * t448 + (-t10 * t341 + t18 * t339 - t9) * t419 + (-t24 * t341 + t25 * t339 - t28) * t420) * t264 + 0.2e1 * ((-qJD(4) * t79 + t31) * t421 + (qJD(4) * t42 - t10 * t335 - t18 * t336 + t303) * t419 + (qJD(4) * t62 - t24 * t335 - t25 * t336 + t302) * t420 + (-t410 / 0.2e1 + t407 / 0.2e1 + t408 / 0.2e1) * t321) * t261 + (-m(5) * t145 - t348) * t250 + (t177 * t250 + t184 * t321 + t197) * m(4); t425 + (-pkin(4) * t28 - t24 * t40 - t25 * t41 - t62 * t80) * m(6) + ((-m(6) * t283 - t260 * t100 - t102 * t263) * qJD(5) - t260 * t46 + t263 * t48 + m(6) * t302) * pkin(10) + t427 * t415 + t429 * t414 + t431 * t400 - t422 * t195 + t230 + t434 * t99 + t435 * t101 + (t1 * t237 + t10 * t435 + t18 * t434 - t2 * t238 + t255 * t9 - t42 * t54) * m(7) + t439 * t385 + ((t416 - Ifges(5,1) / 0.2e1) * t195 + t433 - t458) * t194 + ((t84 + t407) * t383 + t458) * qJD(5) - t358 * t80 + t260 * t452 - pkin(4) * t22 - t54 * t84 - t41 * t100 - t40 * t102 - t79 * t146 - t9 * t299 - t28 * t301 + t302 * mrSges(6,3) + t303 * mrSges(7,3) + t237 * t45 - t238 * t47 + t255 * t21; t277 + (-(-t10 + t17) * t18 + (-t132 * t42 + t1) * pkin(5)) * m(7) + (-t132 * t84 + t45) * pkin(5) + (t10 * t131 + t132 * t18) * mrSges(7,3) + (t131 * t24 + t132 * t25) * mrSges(6,3) + t440 - t17 * t99 - t24 * t100 + t18 * t101 + t25 * t102 - t42 * (mrSges(7,1) * t132 + mrSges(7,2) * t131) - t62 * (mrSges(6,1) * t132 + mrSges(6,2) * t131) + (t131 * t446 - t464) * t402 + t437 * t401 + (t131 * t444 - t132 * t442) * t393 + (-t132 * t443 + t436 + t465) * t404; t132 * t101 - t131 * t99 + 0.2e1 * (t9 / 0.2e1 + t10 * t401 + t18 * t404) * m(7) + t21;];
tauc  = t13(:);

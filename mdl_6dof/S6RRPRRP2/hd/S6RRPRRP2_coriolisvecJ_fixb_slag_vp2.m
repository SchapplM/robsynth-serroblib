% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:06
% EndTime: 2019-03-09 11:43:23
% DurationCPUTime: 8.35s
% Computational Cost: add. (13514->572), mult. (35398->734), div. (0->0), fcn. (26684->8), ass. (0->269)
t406 = Ifges(6,1) + Ifges(7,1);
t404 = Ifges(7,4) + Ifges(6,5);
t405 = Ifges(6,6) - Ifges(7,6);
t248 = sin(qJ(5));
t251 = cos(qJ(5));
t246 = sin(pkin(10));
t247 = cos(pkin(10));
t250 = sin(qJ(2));
t253 = cos(qJ(2));
t220 = -t246 * t250 + t247 * t253;
t209 = t220 * qJD(1);
t319 = qJD(1) * t253;
t320 = qJD(1) * t250;
t210 = -t246 * t319 - t247 * t320;
t249 = sin(qJ(4));
t252 = cos(qJ(4));
t273 = t209 * t249 - t252 * t210;
t315 = qJD(2) + qJD(4);
t144 = t248 * t273 - t251 * t315;
t375 = -t144 / 0.2e1;
t221 = t246 * t253 + t247 * t250;
t211 = t221 * qJD(2);
t196 = qJD(1) * t211;
t212 = t220 * qJD(2);
t197 = qJD(1) * t212;
t120 = qJD(4) * t273 + t252 * t196 + t197 * t249;
t303 = t252 * t209 + t210 * t249;
t119 = qJD(4) * t303 - t196 * t249 + t197 * t252;
t82 = -qJD(5) * t144 + t251 * t119;
t146 = t248 * t315 + t251 * t273;
t83 = qJD(5) * t146 + t248 * t119;
t403 = (-Ifges(6,4) + Ifges(7,5)) * t83 + t406 * t82 + t404 * t120;
t142 = Ifges(6,4) * t144;
t157 = qJD(5) - t303;
t342 = Ifges(7,5) * t144;
t394 = t146 * t406 + t404 * t157 - t142 + t342;
t353 = -qJ(3) - pkin(7);
t231 = t353 * t250;
t225 = qJD(1) * t231;
t232 = t353 * t253;
t226 = qJD(1) * t232;
t325 = t247 * t226;
t170 = -t225 * t246 + t325;
t360 = pkin(8) * t209;
t147 = t170 - t360;
t213 = t246 * t226;
t171 = t247 * t225 + t213;
t359 = pkin(8) * t210;
t148 = t171 + t359;
t240 = pkin(2) * t247 + pkin(3);
t361 = pkin(2) * t246;
t204 = t249 * t240 + t252 * t361;
t391 = -t204 * qJD(4) - t252 * t147 + t148 * t249;
t281 = pkin(5) * t248 - qJ(6) * t251;
t316 = qJD(5) * t251;
t317 = qJD(5) * t248;
t402 = pkin(5) * t317 - qJ(6) * t316 - qJD(6) * t248 - t281 * t303;
t401 = t248 * t404 + t251 * t405;
t340 = Ifges(7,5) * t251;
t343 = Ifges(6,4) * t251;
t400 = t248 * t406 - t340 + t343;
t306 = qJD(2) * t353;
t206 = qJD(3) * t253 + t250 * t306;
t184 = t206 * qJD(1);
t207 = -t250 * qJD(3) + t253 * t306;
t185 = t207 * qJD(1);
t145 = t247 * t184 + t246 * t185;
t127 = -pkin(8) * t196 + t145;
t143 = -t184 * t246 + t185 * t247;
t263 = -pkin(8) * t197 + t143;
t219 = qJD(2) * pkin(2) + t225;
t166 = t247 * t219 + t213;
t136 = qJD(2) * pkin(3) + t166 + t359;
t167 = t246 * t219 - t325;
t140 = t167 + t360;
t88 = t252 * t136 - t249 * t140;
t31 = qJD(4) * t88 + t252 * t127 + t249 * t263;
t243 = pkin(2) * t320;
t238 = qJD(2) * t243;
t173 = pkin(3) * t196 + t238;
t55 = pkin(4) * t120 - pkin(9) * t119 + t173;
t89 = t249 * t136 + t252 * t140;
t86 = pkin(9) * t315 + t89;
t241 = -pkin(2) * t253 - pkin(1);
t321 = qJD(1) * t241;
t227 = qJD(3) + t321;
t172 = -t209 * pkin(3) + t227;
t92 = -pkin(4) * t303 - pkin(9) * t273 + t172;
t6 = t248 * t55 + t251 * t31 + t92 * t316 - t317 * t86;
t38 = t248 * t92 + t251 * t86;
t7 = -qJD(5) * t38 - t248 * t31 + t251 * t55;
t399 = -t7 * t248 + t251 * t6;
t2 = qJ(6) * t120 + qJD(6) * t157 + t6;
t37 = -t248 * t86 + t251 * t92;
t388 = qJD(6) - t37;
t26 = -pkin(5) * t157 + t388;
t27 = qJ(6) * t157 + t38;
t4 = -pkin(5) * t120 - t7;
t398 = t2 * t251 + t248 * t4 + t26 * t316 - t27 * t317;
t290 = -Ifges(6,2) * t248 + t343;
t295 = mrSges(7,1) * t248 - mrSges(7,3) * t251;
t297 = mrSges(6,1) * t248 + mrSges(6,2) * t251;
t85 = -pkin(4) * t315 - t88;
t44 = t144 * pkin(5) - t146 * qJ(6) + t85;
t397 = t290 * t375 + t85 * t297 + t44 * t295 - (t248 * t38 + t251 * t37) * mrSges(6,3);
t374 = t144 / 0.2e1;
t373 = -t146 / 0.2e1;
t370 = -t157 / 0.2e1;
t339 = Ifges(5,2) * t303;
t396 = t339 / 0.2e1;
t395 = pkin(5) * t273;
t393 = -t89 + t402;
t392 = t402 - t391;
t335 = -mrSges(5,1) * t315 + mrSges(6,1) * t144 + mrSges(6,2) * t146 + mrSges(5,3) * t273;
t390 = qJ(6) * t273;
t169 = t220 * t249 + t221 * t252;
t186 = -t220 * pkin(3) + t241;
t272 = t252 * t220 - t221 * t249;
t109 = -pkin(4) * t272 - t169 * pkin(9) + t186;
t175 = t247 * t231 + t232 * t246;
t154 = -pkin(8) * t221 + t175;
t176 = t246 * t231 - t247 * t232;
t155 = pkin(8) * t220 + t176;
t111 = t154 * t249 + t155 * t252;
t387 = t248 * t109 + t251 * t111;
t386 = t252 * t154 - t155 * t249;
t203 = t240 * t252 - t249 * t361;
t156 = Ifges(5,4) * t303;
t286 = Ifges(6,5) * t251 - Ifges(6,6) * t248;
t264 = t157 * t286;
t288 = Ifges(7,4) * t251 + Ifges(7,6) * t248;
t265 = t157 * t288;
t341 = Ifges(7,5) * t248;
t292 = Ifges(7,1) * t251 + t341;
t266 = t146 * t292;
t344 = Ifges(6,4) * t248;
t294 = Ifges(6,1) * t251 - t344;
t267 = t146 * t294;
t284 = Ifges(7,3) * t248 + t340;
t268 = t144 * t284;
t363 = -t251 / 0.2e1;
t364 = t248 / 0.2e1;
t365 = -t248 / 0.2e1;
t348 = Ifges(5,1) * t273;
t383 = t156 / 0.2e1 + t348 / 0.2e1;
t141 = Ifges(7,5) * t146;
t71 = t157 * Ifges(7,6) + t144 * Ifges(7,3) + t141;
t345 = Ifges(6,4) * t146;
t74 = -t144 * Ifges(6,2) + t157 * Ifges(6,6) + t345;
t254 = (t248 * t27 - t251 * t26) * mrSges(7,2) + t88 * mrSges(5,3) - Ifges(5,5) * t315 - t268 / 0.2e1 - t267 / 0.2e1 - t266 / 0.2e1 - t265 / 0.2e1 - t264 / 0.2e1 - t172 * mrSges(5,2) + t71 * t365 + t74 * t364 + t394 * t363 - t383 - t397;
t385 = t254 - t156 / 0.2e1;
t349 = mrSges(6,3) * t146;
t105 = mrSges(6,1) * t157 - t349;
t106 = -mrSges(7,1) * t157 + mrSges(7,2) * t146;
t322 = t105 - t106;
t103 = -mrSges(7,2) * t144 + mrSges(7,3) * t157;
t350 = mrSges(6,3) * t144;
t104 = -mrSges(6,2) * t157 - t350;
t323 = t103 + t104;
t41 = mrSges(6,1) * t120 - mrSges(6,3) * t82;
t42 = -t120 * mrSges(7,1) + t82 * mrSges(7,2);
t351 = -t41 + t42;
t40 = -mrSges(7,2) * t83 + mrSges(7,3) * t120;
t43 = -mrSges(6,2) * t120 - mrSges(6,3) * t83;
t352 = t40 + t43;
t384 = m(6) * (-t316 * t37 - t317 * t38 + t399) + m(7) * t398 + t352 * t251 + t351 * t248 + (-t248 * t323 - t251 * t322) * qJD(5);
t152 = -t206 * t246 + t247 * t207;
t133 = -pkin(8) * t212 + t152;
t153 = t247 * t206 + t246 * t207;
t134 = -pkin(8) * t211 + t153;
t50 = qJD(4) * t386 + t133 * t249 + t134 * t252;
t124 = qJD(4) * t272 - t211 * t249 + t212 * t252;
t125 = qJD(4) * t169 + t252 * t211 + t212 * t249;
t318 = qJD(2) * t250;
t180 = pkin(2) * t318 + pkin(3) * t211;
t61 = pkin(4) * t125 - pkin(9) * t124 + t180;
t13 = -qJD(5) * t387 - t248 * t50 + t251 * t61;
t123 = pkin(4) * t273 - pkin(9) * t303;
t312 = -Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1;
t313 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t314 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t382 = t313 * t144 + t314 * t146 - t312 * t157 + t172 * mrSges(5,1) + t37 * mrSges(6,1) - t26 * mrSges(7,1) - t38 * mrSges(6,2) - t89 * mrSges(5,3) + t27 * mrSges(7,3) - Ifges(5,4) * t273 - Ifges(5,6) * t315 - Ifges(6,6) * t374 - Ifges(7,6) * t375 - t396 - t404 * t373 - (Ifges(7,2) + Ifges(6,3)) * t370;
t381 = t82 / 0.2e1;
t380 = -t83 / 0.2e1;
t379 = t83 / 0.2e1;
t378 = pkin(1) * mrSges(3,1);
t377 = pkin(1) * mrSges(3,2);
t376 = t120 / 0.2e1;
t372 = t146 / 0.2e1;
t368 = -t210 / 0.2e1;
t367 = -t211 / 0.2e1;
t366 = t212 / 0.2e1;
t362 = t251 / 0.2e1;
t95 = t147 * t249 + t148 * t252;
t179 = -pkin(3) * t210 + t243;
t96 = t123 + t179;
t49 = t248 * t96 + t251 * t95;
t347 = Ifges(3,4) * t250;
t32 = qJD(4) * t89 + t127 * t249 - t252 * t263;
t338 = t386 * t32;
t337 = t210 * Ifges(4,4);
t54 = t248 * t123 + t251 * t88;
t334 = Ifges(3,5) * qJD(2);
t333 = Ifges(3,6) * qJD(2);
t332 = qJD(2) * mrSges(3,1);
t331 = qJD(2) * mrSges(3,2);
t191 = t203 * qJD(4);
t329 = t191 * t248;
t328 = t191 * t251;
t309 = t334 / 0.2e1;
t308 = -t333 / 0.2e1;
t305 = t196 * mrSges(4,1) + t197 * mrSges(4,2);
t304 = t120 * mrSges(5,1) + t119 * mrSges(5,2);
t300 = -t2 * t248 + t251 * t4;
t299 = -t248 * t6 - t251 * t7;
t298 = mrSges(6,1) * t251 - mrSges(6,2) * t248;
t296 = mrSges(7,1) * t251 + mrSges(7,3) * t248;
t289 = Ifges(6,2) * t251 + t344;
t283 = -Ifges(7,3) * t251 + t341;
t282 = pkin(5) * t251 + qJ(6) * t248;
t279 = t248 * t26 + t251 * t27;
t276 = t248 * t37 - t251 * t38;
t48 = -t248 * t95 + t251 * t96;
t53 = t123 * t251 - t248 * t88;
t58 = t109 * t251 - t111 * t248;
t228 = -pkin(4) - t282;
t12 = t109 * t316 - t111 * t317 + t248 * t61 + t251 * t50;
t262 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t51 = qJD(4) * t111 - t252 * t133 + t134 * t249;
t11 = pkin(5) * t83 - qJ(6) * t82 - qJD(6) * t146 + t32;
t19 = t82 * Ifges(7,5) + t120 * Ifges(7,6) + t83 * Ifges(7,3);
t20 = t82 * Ifges(6,4) - t83 * Ifges(6,2) + t120 * Ifges(6,6);
t256 = -t31 * mrSges(5,2) + Ifges(5,5) * t119 - Ifges(5,6) * t120 - t11 * t296 + t19 * t363 + t20 * t362 + t283 * t379 + t289 * t380 + t400 * t381 + t401 * t376 + t403 * t364 + (-mrSges(5,1) - t298) * t32 + (-t74 / 0.2e1 + t71 / 0.2e1) * t317 + t394 * t316 / 0.2e1 + t399 * mrSges(6,3) + t397 * qJD(5) + t398 * mrSges(7,2) + (t268 + t267 + t266 + t265 + t264) * qJD(5) / 0.2e1;
t242 = Ifges(3,4) * t319;
t230 = mrSges(3,3) * t319 - t331;
t229 = -mrSges(3,3) * t320 + t332;
t218 = Ifges(3,1) * t320 + t242 + t334;
t217 = t333 + (Ifges(3,2) * t253 + t347) * qJD(1);
t202 = Ifges(4,4) * t209;
t200 = -pkin(4) - t203;
t183 = qJD(2) * mrSges(4,1) + t210 * mrSges(4,3);
t182 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t209;
t174 = -t203 + t228;
t165 = -mrSges(4,1) * t209 - mrSges(4,2) * t210;
t159 = -t210 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t202;
t158 = t209 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t337;
t149 = -mrSges(5,2) * t315 + mrSges(5,3) * t303;
t122 = -mrSges(5,1) * t303 + mrSges(5,2) * t273;
t117 = Ifges(7,2) * t120;
t115 = Ifges(6,3) * t120;
t98 = mrSges(7,1) * t144 - mrSges(7,3) * t146;
t97 = pkin(5) * t146 + qJ(6) * t144;
t81 = Ifges(7,4) * t82;
t80 = Ifges(6,5) * t82;
t79 = Ifges(6,6) * t83;
t78 = Ifges(7,6) * t83;
t62 = t169 * t281 - t386;
t46 = pkin(5) * t272 - t58;
t45 = -qJ(6) * t272 + t387;
t36 = -t53 - t395;
t35 = t54 + t390;
t34 = -t48 - t395;
t33 = t49 + t390;
t24 = mrSges(6,1) * t83 + mrSges(6,2) * t82;
t23 = mrSges(7,1) * t83 - mrSges(7,3) * t82;
t14 = t281 * t124 + (qJD(5) * t282 - qJD(6) * t251) * t169 + t51;
t9 = -pkin(5) * t125 - t13;
t8 = qJ(6) * t125 - qJD(6) * t272 + t12;
t1 = [t227 * (mrSges(4,1) * t211 + mrSges(4,2) * t212) + t8 * t103 + t12 * t104 + t13 * t105 + t9 * t106 + t14 * t98 + m(4) * (t143 * t175 + t145 * t176 + t152 * t166 + t153 * t167) + t62 * t23 + t58 * t41 + t45 * t40 + t46 * t42 + (-t254 + t383) * t124 + (-t339 / 0.2e1 + t382) * t125 + (-t143 * t221 + t145 * t220 - t166 * t212 - t167 * t211 - t175 * t197 - t176 * t196) * mrSges(4,3) + (-t221 * t196 + t220 * t197 + t209 * t366 - t211 * t368) * Ifges(4,4) + (-t220 * t196 + t209 * t367) * Ifges(4,2) + m(6) * (t12 * t38 + t13 * t37 + t387 * t6 + t51 * t85 + t58 * t7 - t338) + m(5) * (t111 * t31 + t172 * t180 + t173 * t186 + t50 * t89 - t51 * t88 - t338) + t387 * t43 - (t80 / 0.2e1 - t79 / 0.2e1 + t115 / 0.2e1 + t81 / 0.2e1 + t117 / 0.2e1 + t78 / 0.2e1 - Ifges(5,4) * t119 + t173 * mrSges(5,1) - t31 * mrSges(5,3) + t313 * t83 + t314 * t82 + (Ifges(5,2) - t312) * t120 + t262) * t272 + m(7) * (t11 * t62 + t14 * t44 + t2 * t45 + t26 * t9 + t27 * t8 + t4 * t46) + t159 * t366 + t158 * t367 - t386 * t24 + (-t111 * t120 - t119 * t386) * mrSges(5,3) + (t11 * t295 + t284 * t379 + t290 * t380 + Ifges(5,1) * t119 - Ifges(5,4) * t120 + t19 * t364 + t173 * mrSges(5,2) + (mrSges(5,3) + t297) * t32 + t299 * mrSges(6,3) + t300 * mrSges(7,2) + (-mrSges(7,2) * t279 + mrSges(6,3) * t276 + t283 * t375 + t289 * t374 + t296 * t44 + t298 * t85 + t363 * t74 + t370 * t401 + t373 * t400) * qJD(5) + (t292 + t294) * t381 + (t288 + t286) * t376 + (qJD(5) * t394 + t20) * t365 + (qJD(5) * t71 + t403) * t362) * t169 + t50 * t149 + t180 * t122 + t153 * t182 + t152 * t183 + t186 * t304 + t241 * t305 + t335 * t51 + (t221 * t197 + t212 * t368) * Ifges(4,1) + (Ifges(4,5) * t366 + Ifges(4,6) * t367 + (-pkin(7) * t229 + t218 / 0.2e1 + t309 + (-0.2e1 * t377 + 0.3e1 / 0.2e1 * Ifges(3,4) * t253) * qJD(1)) * t253) * qJD(2) + (-pkin(7) * t230 - t217 / 0.2e1 + t308 + (-0.2e1 * t378 - 0.3e1 / 0.2e1 * t347 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t253) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t220 + mrSges(4,2) * t221) + t165 + m(4) * (t227 + t321)) * pkin(2)) * t318; -(Ifges(4,2) * t210 + t159 + t202) * t209 / 0.2e1 - t33 * t103 - t49 * t104 - t48 * t105 - t34 * t106 + m(4) * (t143 * t247 + t145 * t246) * pkin(2) + t256 + (-t119 * t203 - t120 * t204) * mrSges(5,3) - m(4) * (t166 * t170 + t167 * t171) + (t11 * t174 + t392 * t44 + (t328 - t33) * t27 + (t329 - t34) * t26) * m(7) + t392 * t98 + (-t172 * t179 - t203 * t32 + t204 * t31 + (-t95 + t191) * t89 + t391 * t88) * m(5) - t335 * t391 + (t200 * t32 - t391 * t85 + (t328 - t49) * t38 + (-t329 - t48) * t37) * m(6) + (-t348 / 0.2e1 + t385) * t303 + t384 * (pkin(9) + t204) + (t166 * t209 - t167 * t210 + (-t196 * t246 - t197 * t247) * pkin(2)) * mrSges(4,3) + t158 * t368 + (-t382 + t396) * t273 + t143 * mrSges(4,1) - t145 * mrSges(4,2) - t95 * t149 + t174 * t23 - t179 * t122 - t171 * t182 - t170 * t183 - Ifges(4,6) * t196 + Ifges(4,5) * t197 + t200 * t24 - qJD(2) * (Ifges(4,5) * t209 + Ifges(4,6) * t210) / 0.2e1 - t227 * (-mrSges(4,1) * t210 + mrSges(4,2) * t209) + (-t248 * t322 + t251 * t323 + t149) * t191 + t210 * (Ifges(4,1) * t209 + t337) / 0.2e1 + ((t309 - t242 / 0.2e1 - t218 / 0.2e1 + qJD(1) * t377 + (t229 - t332) * pkin(7)) * t253 + (t308 + t217 / 0.2e1 + (t378 + t347 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t253) * qJD(1) + (t230 + t331) * pkin(7) + (-m(4) * t227 - t165) * pkin(2)) * t250) * qJD(1); -t303 * t149 - t209 * t182 - t210 * t183 + (-t98 - t335) * t273 + (t157 * t323 - t351) * t251 + (-t157 * t322 + t352) * t248 + t304 + t305 + (t157 * t279 - t273 * t44 - t300) * m(7) + (-t157 * t276 - t273 * t85 - t299) * m(6) + (t273 * t88 - t303 * t89 + t173) * m(5) + (-t166 * t210 - t167 * t209 + t238) * m(4); ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t273 + t385) * t303 - t35 * t103 - t54 * t104 - t53 * t105 - t36 * t106 + t256 - t335 * t89 - pkin(4) * t24 + t384 * pkin(9) + t393 * t98 - t88 * t149 - t382 * t273 + t228 * t23 + (t11 * t228 - t26 * t36 - t27 * t35 + t393 * t44) * m(7) + (-pkin(4) * t32 - t37 * t53 - t38 * t54 - t85 * t89) * m(6); (t144 * t26 + t146 * t27) * mrSges(7,2) + qJD(6) * t103 - t97 * t98 + t117 + t115 + t81 + t80 - t79 + t78 - pkin(5) * t42 + qJ(6) * t40 + (Ifges(7,3) * t146 - t342) * t375 + t74 * t372 - t44 * (mrSges(7,1) * t146 + mrSges(7,3) * t144) - t85 * (mrSges(6,1) * t146 - mrSges(6,2) * t144) + t262 + (t322 + t349) * t38 + (-t323 - t350) * t37 + (-t404 * t144 - t146 * t405) * t370 + (-pkin(5) * t4 + qJ(6) * t2 - t26 * t38 + t27 * t388 - t44 * t97) * m(7) + (-Ifges(6,2) * t146 - t142 + t394) * t374 + (-t144 * t406 + t141 - t345 + t71) * t373; -t157 * t103 + t146 * t98 + 0.2e1 * (t4 / 0.2e1 + t44 * t372 + t27 * t370) * m(7) + t42;];
tauc  = t1(:);

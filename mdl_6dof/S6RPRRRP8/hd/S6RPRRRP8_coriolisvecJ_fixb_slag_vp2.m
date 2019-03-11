% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:57
% EndTime: 2019-03-09 06:23:13
% DurationCPUTime: 9.25s
% Computational Cost: add. (7643->518), mult. (16444->686), div. (0->0), fcn. (10473->6), ass. (0->242)
t391 = Ifges(6,1) + Ifges(7,1);
t405 = Ifges(6,4) - Ifges(7,5);
t390 = Ifges(7,4) + Ifges(6,5);
t389 = Ifges(7,2) + Ifges(6,3);
t388 = Ifges(6,6) - Ifges(7,6);
t220 = sin(qJ(4));
t349 = cos(qJ(4));
t350 = cos(qJ(3));
t262 = t349 * t350;
t298 = qJD(3) + qJD(4);
t227 = t298 * t262;
t221 = sin(qJ(3));
t299 = qJD(1) * qJD(3);
t272 = t221 * t299;
t305 = qJD(1) * t221;
t279 = t220 * t305;
t144 = qJD(1) * t227 - qJD(4) * t279 - t220 * t272;
t282 = t349 * t221;
t285 = t220 * t350;
t235 = -t282 - t285;
t178 = t235 * qJD(1);
t174 = qJD(5) - t178;
t223 = -pkin(1) - pkin(7);
t197 = qJD(1) * t223 + qJD(2);
t169 = -pkin(8) * t305 + t197 * t221;
t165 = t349 * t169;
t187 = t350 * t197;
t275 = qJD(1) * t350;
t170 = -pkin(8) * t275 + t187;
t166 = qJD(3) * pkin(3) + t170;
t117 = t220 * t166 + t165;
t106 = pkin(9) * t298 + t117;
t179 = qJD(1) * t262 - t279;
t194 = pkin(3) * t305 + qJD(1) * qJ(2);
t118 = -pkin(4) * t178 - pkin(9) * t179 + t194;
t219 = sin(qJ(5));
t222 = cos(qJ(5));
t300 = qJD(5) * t222;
t301 = qJD(5) * t219;
t164 = t220 * t169;
t116 = t349 * t166 - t164;
t167 = t170 * qJD(3);
t304 = qJD(3) * t221;
t238 = (pkin(8) * qJD(1) - t197) * t304;
t43 = qJD(4) * t116 + t349 * t167 + t220 * t238;
t143 = t298 * t178;
t217 = qJD(1) * qJD(2);
t274 = qJD(3) * t350;
t260 = qJD(1) * t274;
t188 = pkin(3) * t260 + t217;
t56 = pkin(4) * t144 - pkin(9) * t143 + t188;
t6 = -t106 * t301 + t118 * t300 + t219 * t56 + t222 * t43;
t2 = qJ(6) * t144 + qJD(6) * t174 + t6;
t40 = t106 * t222 + t118 * t219;
t7 = -qJD(5) * t40 - t219 * t43 + t222 * t56;
t4 = -pkin(5) * t144 - t7;
t257 = t2 * t222 + t219 * t4;
t39 = -t106 * t219 + t118 * t222;
t381 = qJD(6) - t39;
t26 = -pkin(5) * t174 + t381;
t404 = t26 * t300 + t257;
t158 = t219 * t179 - t222 * t298;
t302 = qJD(5) * t158;
t91 = t222 * t143 - t302;
t159 = t222 * t179 + t219 * t298;
t92 = qJD(5) * t159 + t219 * t143;
t403 = t390 * t144 + t391 * t91 - t405 * t92;
t402 = -t158 * t388 + t159 * t390 + t389 * t174;
t157 = Ifges(6,4) * t158;
t332 = Ifges(7,5) * t158;
t387 = t159 * t391 + t390 * t174 - t157 + t332;
t365 = -t92 / 0.2e1;
t401 = Ifges(6,2) * t365;
t324 = t179 * mrSges(5,3);
t306 = mrSges(5,1) * t298 - mrSges(6,1) * t158 - mrSges(6,2) * t159 - t324;
t315 = t178 * t222;
t316 = t178 * t219;
t400 = -qJD(6) * t219 + (-t300 + t315) * qJ(6) + (t301 - t316) * pkin(5);
t399 = -t219 * t388 + t222 * t390;
t331 = Ifges(7,5) * t219;
t334 = Ifges(6,4) * t219;
t398 = t222 * t391 + t331 - t334;
t120 = t170 * t220 + t165;
t303 = qJD(4) * t220;
t397 = pkin(3) * t303 - t120;
t363 = t144 / 0.2e1;
t366 = t91 / 0.2e1;
t396 = Ifges(6,4) * t366 + Ifges(6,6) * t363 - t91 * Ifges(7,5) / 0.2e1 - t144 * Ifges(7,6) / 0.2e1 + Ifges(7,3) * t365 + t401;
t29 = -mrSges(7,2) * t92 + mrSges(7,3) * t144;
t30 = mrSges(6,1) * t144 - mrSges(6,3) * t91;
t31 = -t144 * mrSges(7,1) + t91 * mrSges(7,2);
t32 = -mrSges(6,2) * t144 - mrSges(6,3) * t92;
t395 = (-t30 + t31) * t219 + (t29 + t32) * t222;
t394 = (m(4) + m(3)) * qJ(2) + mrSges(3,3);
t273 = qJD(4) * t349;
t152 = -qJD(3) * t282 - qJD(4) * t285 - t220 * t274 - t221 * t273;
t392 = t152 / 0.2e1;
t355 = -t178 / 0.2e1;
t340 = pkin(8) - t223;
t255 = mrSges(7,1) * t219 - mrSges(7,3) * t222;
t105 = -pkin(4) * t298 - t116;
t35 = t158 * pkin(5) - t159 * qJ(6) + t105;
t386 = t255 * t35;
t385 = t397 + t400;
t384 = -t117 + t400;
t24 = mrSges(6,1) * t92 + mrSges(6,2) * t91;
t328 = t143 * mrSges(5,3);
t383 = t328 + t24;
t256 = mrSges(6,1) * t219 + mrSges(6,2) * t222;
t382 = t105 * t256;
t184 = t220 * t221 - t262;
t208 = t221 * pkin(3) + qJ(2);
t148 = -pkin(4) * t235 + pkin(9) * t184 + t208;
t189 = t340 * t221;
t190 = t340 * t350;
t155 = -t189 * t349 - t220 * t190;
t380 = t219 * t148 + t222 * t155;
t379 = t220 * t189 - t349 * t190;
t348 = pkin(3) * t220;
t210 = pkin(9) + t348;
t265 = pkin(3) * t273;
t378 = -t210 * t301 + t222 * t265;
t377 = -t210 * t300 - t219 * t265;
t375 = t389 * t144 - t388 * t92 + t390 * t91;
t343 = t219 * t7;
t374 = -t39 * t300 - t40 * t301 - t343;
t153 = -t220 * t304 - t221 * t303 + t227;
t44 = qJD(4) * t117 + t220 * t167 - t349 * t238;
t323 = t184 * t44;
t373 = t116 * t152 + t117 * t153 + t323;
t198 = pkin(3) * t274 + qJD(2);
t82 = pkin(4) * t153 - pkin(9) * t152 + t198;
t182 = t340 * t304;
t183 = t190 * qJD(3);
t93 = qJD(4) * t379 + t220 * t182 - t349 * t183;
t13 = -qJD(5) * t380 - t219 * t93 + t222 * t82;
t370 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t369 = qJD(1) ^ 2;
t364 = t92 / 0.2e1;
t361 = -t158 / 0.2e1;
t360 = t158 / 0.2e1;
t359 = -t159 / 0.2e1;
t358 = t159 / 0.2e1;
t357 = -t174 / 0.2e1;
t353 = t179 / 0.2e1;
t347 = pkin(5) * t179;
t9 = t92 * pkin(5) - t91 * qJ(6) - t159 * qJD(6) + t44;
t346 = t184 * t9;
t342 = t222 * t6;
t339 = mrSges(6,3) * t158;
t338 = mrSges(6,3) * t159;
t337 = Ifges(4,4) * t221;
t336 = Ifges(5,4) * t179;
t335 = Ifges(6,4) * t159;
t333 = Ifges(6,4) * t222;
t330 = Ifges(7,5) * t222;
t329 = pkin(3) * qJD(4);
t327 = t144 * mrSges(5,3);
t326 = t379 * t44;
t325 = t178 * mrSges(5,3);
t320 = t152 * t219;
t319 = t152 * t222;
t318 = t153 * t219;
t317 = t153 * t222;
t313 = t184 * t222;
t193 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t275;
t309 = t221 * t193;
t147 = t179 * pkin(4) - t178 * pkin(9);
t55 = t222 * t116 + t219 * t147;
t109 = -mrSges(7,2) * t158 + mrSges(7,3) * t174;
t110 = -mrSges(6,2) * t174 - t339;
t308 = t109 + t110;
t111 = mrSges(6,1) * t174 - t338;
t112 = -mrSges(7,1) * t174 + mrSges(7,2) * t159;
t307 = -t111 + t112;
t121 = t170 * t349 - t164;
t264 = pkin(3) * t275;
t127 = t264 + t147;
t52 = t222 * t121 + t219 * t127;
t295 = t349 * pkin(3);
t292 = Ifges(4,4) * t350;
t27 = qJ(6) * t174 + t40;
t291 = t27 * t301;
t78 = -t158 * Ifges(6,2) + t174 * Ifges(6,6) + t335;
t287 = -t219 * t78 / 0.2e1;
t286 = t219 * t349;
t284 = t222 * t349;
t268 = t301 / 0.2e1;
t267 = t300 / 0.2e1;
t252 = -Ifges(6,2) * t219 + t333;
t249 = Ifges(7,3) * t219 + t330;
t248 = -t222 * pkin(5) - t219 * qJ(6);
t247 = pkin(5) * t219 - qJ(6) * t222;
t246 = t219 * t27 - t222 * t26;
t245 = t219 * t40 + t222 * t39;
t54 = -t116 * t219 + t147 * t222;
t51 = -t121 * t219 + t127 * t222;
t83 = t148 * t222 - t155 * t219;
t191 = -pkin(4) + t248;
t239 = -Ifges(4,5) * t221 - Ifges(4,6) * t350;
t237 = t184 * t300 - t320;
t236 = t184 * t301 + t319;
t12 = t148 * t300 - t155 * t301 + t219 * t82 + t222 * t93;
t234 = qJ(2) * (mrSges(4,1) * t350 - mrSges(4,2) * t221);
t233 = t221 * (-Ifges(4,2) * t350 - t337);
t232 = (Ifges(4,1) * t350 - t337) * qJD(1);
t231 = (-Ifges(4,2) * t221 + t292) * qJD(1);
t186 = (t221 * mrSges(4,1) + mrSges(4,2) * t350) * qJD(1);
t230 = -t219 * t308 + t222 * t307;
t229 = (-Ifges(4,1) * t221 - t292) * t350;
t228 = -qJD(5) * t245 - t343;
t94 = qJD(4) * t155 - t349 * t182 - t220 * t183;
t226 = t230 * qJD(5) + m(6) * (t342 + t374) + m(7) * (-t291 + t404) + t395;
t128 = Ifges(5,2) * t178 + t298 * Ifges(5,6) + t336;
t172 = Ifges(5,4) * t178;
t129 = Ifges(5,1) * t179 + t298 * Ifges(5,5) + t172;
t156 = Ifges(7,5) * t159;
t75 = t174 * Ifges(7,6) + t158 * Ifges(7,3) + t156;
t225 = (-t26 * t315 + t404) * mrSges(7,2) - t178 * t386 - t178 * t382 + (t403 / 0.2e1 + mrSges(6,2) * t44 - mrSges(7,3) * t9 + t390 * t363 + t391 * t366) * t219 + (Ifges(6,6) * t179 + t178 * t252) * t360 + (Ifges(7,6) * t179 + t178 * t249) * t361 + (-Ifges(5,2) * t179 + t129 + t172) * t355 + (t382 + t386 + t287) * qJD(5) + mrSges(6,3) * t342 + t128 * t353 + t116 * t325 + (-t315 / 0.2e1 + t267) * t387 + (-t252 / 0.2e1 + t249 / 0.2e1) * t302 - (Ifges(5,1) * t178 - t336 + t402) * t179 / 0.2e1 + (-mrSges(6,1) * t44 - mrSges(7,1) * t9 - Ifges(7,3) * t364 + t363 * t388 + t396 + t401) * t222 + (t178 * t399 + t179 * t389) * t357 + (t178 * t398 + t179 * t390) * t359 + t331 * t364 + t334 * t365 - t298 * (Ifges(5,5) * t178 - Ifges(5,6) * t179) / 0.2e1 + (-t330 + t333) * t366 + (t40 * mrSges(6,3) + t78 / 0.2e1 - t75 / 0.2e1) * t316 + t40 * mrSges(6,2) * t179 + t26 * mrSges(7,1) * t179 - t39 * (mrSges(6,1) * t179 - mrSges(6,3) * t315) - t194 * (mrSges(5,1) * t179 + mrSges(5,2) * t178) + (t159 * t398 + t174 * t399) * qJD(5) / 0.2e1 + Ifges(5,5) * t143 - Ifges(5,6) * t144 - t43 * mrSges(5,2) - t44 * mrSges(5,1) + t75 * t268;
t211 = -t295 - pkin(4);
t192 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t305;
t180 = -t295 + t191;
t177 = Ifges(4,5) * qJD(3) + t232;
t176 = Ifges(4,6) * qJD(3) + t231;
t171 = t179 * qJ(6);
t161 = -mrSges(5,2) * t298 + t325;
t146 = -mrSges(5,1) * t178 + mrSges(5,2) * t179;
t138 = -mrSges(7,2) * t316 + mrSges(7,3) * t179;
t101 = mrSges(7,1) * t158 - mrSges(7,3) * t159;
t100 = pkin(5) * t159 + qJ(6) * t158;
t97 = -t184 * t247 - t379;
t63 = pkin(5) * t235 - t83;
t62 = -qJ(6) * t235 + t380;
t38 = -t54 - t347;
t37 = t171 + t55;
t34 = -t51 - t347;
t33 = t171 + t52;
t23 = mrSges(7,1) * t92 - mrSges(7,3) * t91;
t14 = t247 * t152 + (qJD(5) * t248 + qJD(6) * t222) * t184 + t94;
t11 = -pkin(5) * t153 - t13;
t10 = qJ(6) * t153 - qJD(6) * t235 + t12;
t1 = [(-m(5) * t116 + m(6) * t105 - t306) * t94 + (-t153 * t353 + t178 * t392) * Ifges(5,4) + t153 * t355 * Ifges(5,2) - t373 * mrSges(5,3) + (-mrSges(7,2) * t4 + mrSges(6,3) * t7) * t313 + (0.2e1 * t234 - t233 + t229) * t299 + (t390 * t153 + t391 * t236 + t237 * t405) * t358 + (t192 * t274 - t193 * t304) * t223 + (Ifges(7,5) * t236 + Ifges(7,6) * t153 - Ifges(7,3) * t237) * t360 + (Ifges(6,4) * t236 + Ifges(6,2) * t237 + Ifges(6,6) * t153) * t361 + t26 * (-mrSges(7,1) * t153 + mrSges(7,2) * t236) + t39 * (mrSges(6,1) * t153 - mrSges(6,3) * t236) + t40 * (-mrSges(6,2) * t153 + mrSges(6,3) * t237) + t27 * (mrSges(7,2) * t237 + mrSges(7,3) * t153) + t105 * (-mrSges(6,1) * t237 + mrSges(6,2) * t236) + t35 * (-mrSges(7,1) * t237 - mrSges(7,3) * t236) + t129 * t392 - t255 * t346 - t155 * t327 - t256 * t323 + t380 * t32 + (-t375 / 0.2e1 + t43 * mrSges(5,3) - t188 * mrSges(5,1) - Ifges(6,6) * t365 - Ifges(7,6) * t364 - t363 * t389 - t366 * t390 - t370 - Ifges(5,2) * t144 + Ifges(5,4) * t143) * t235 + m(5) * (t117 * t93 + t155 * t43 + t188 * t208 + t194 * t198 - t326) + m(6) * (t12 * t40 + t13 * t39 + t380 * t6 + t7 * t83 - t326) + t402 * t153 / 0.2e1 - (qJD(5) * t75 + t403) * t313 / 0.2e1 + t298 * (Ifges(5,5) * t152 - Ifges(5,6) * t153) / 0.2e1 + t152 * t287 - (t231 + t176) * t274 / 0.2e1 - (t232 + t177) * t304 / 0.2e1 - t383 * t379 + t75 * t320 / 0.2e1 + t208 * (mrSges(5,1) * t144 + mrSges(5,2) * t143) + t194 * (mrSges(5,1) * t153 + mrSges(5,2) * t152) + t198 * t146 + 0.2e1 * qJD(2) * t186 + Ifges(5,1) * t152 * t353 + t93 * t161 - t153 * t128 / 0.2e1 + t10 * t109 + t12 * t110 + t13 * t111 + t11 * t112 + t14 * t101 + t97 * t23 + t83 * t30 + qJD(3) ^ 2 * t239 / 0.2e1 + t62 * t29 + t63 * t31 + t387 * t319 / 0.2e1 + (t153 * t389 + t236 * t390 + t237 * t388) * t174 / 0.2e1 + m(7) * (t10 * t27 + t11 * t26 + t14 * t35 + t2 * t62 + t4 * t63 + t9 * t97) + 0.2e1 * t394 * t217 + (-t188 * mrSges(5,2) - t249 * t364 - t252 * t365 + t267 * t78 - t398 * t366 - t399 * t363 - Ifges(5,1) * t143 + t387 * t268 + (mrSges(7,2) * t2 + mrSges(6,3) * t6 + t396) * t219 + Ifges(5,4) * t144) * t184; (t192 * t350 - t309) * qJD(3) - t394 * t369 + (t23 + t383) * t184 + (-t101 + t306) * t152 + (t219 * t307 + t222 * t308 + t161) * t153 + m(5) * t373 + m(7) * (-t152 * t35 + t26 * t318 + t27 * t317 + t346) + m(6) * (-t105 * t152 + t317 * t40 - t318 * t39 + t323) - (m(5) * t43 + t226 - t327) * t235 + (-m(5) * t194 - m(6) * t245 - m(7) * t246 - t146 - t186 + t230) * qJD(1); t385 * t101 + (-t34 - t377) * t112 + (-t51 + t377) * t111 + (-t52 + t378) * t110 + (-t33 + t378) * t109 + t374 * mrSges(6,3) + (t233 / 0.2e1 - t234 - t229 / 0.2e1) * t369 + t225 + t177 * t305 / 0.2e1 - mrSges(7,2) * t291 + (t265 - t121) * t161 - t192 * t187 + t117 * t324 + (-mrSges(4,1) * t304 - mrSges(4,2) * t274 + t309) * t197 - t295 * t328 - t327 * t348 + t176 * t275 / 0.2e1 - t239 * t299 / 0.2e1 - Ifges(4,6) * t260 + ((-t349 * t44 + t220 * t43 + (-t116 * t220 + t117 * t349) * qJD(4)) * pkin(3) + t116 * t120 - t117 * t121 - t194 * t264) * m(5) + (-t26 * t34 - t27 * t33 + t9 * t180 + (t26 * t286 + t27 * t284) * t329 + t385 * t35) * m(7) + (-t105 * t120 - t39 * t51 - t40 * t52 + t44 * t211 + (t105 * t220 + t284 * t40 - t286 * t39) * t329) * m(6) + t211 * t24 + t180 * t23 - Ifges(4,5) * t272 - t27 * t138 - t146 * t264 + ((-qJD(5) * t246 + t257) * m(7) + (t228 + t342) * m(6) + t395) * t210 - t306 * t397; t225 + (-mrSges(7,2) * t301 - t138) * t27 + (t306 + t324) * t117 + t384 * t101 + t228 * mrSges(6,3) + t226 * pkin(9) + t191 * t23 - t116 * t161 - t37 * t109 - t55 * t110 - t54 * t111 - t38 * t112 - pkin(4) * t24 + (t191 * t9 - t26 * t38 - t27 * t37 + t35 * t384) * m(7) + (-pkin(4) * t44 - t105 * t117 - t39 * t54 - t40 * t55) * m(6); t370 + (Ifges(7,3) * t159 - t332) * t361 + (-pkin(5) * t4 + qJ(6) * t2 - t100 * t35 - t26 * t40 + t27 * t381) * m(7) + (-t158 * t390 - t159 * t388) * t357 + (-t307 + t338) * t40 + (-t308 - t339) * t39 + (-t158 * t391 + t156 - t335 + t75) * t359 + t78 * t358 + (-Ifges(6,2) * t159 - t157 + t387) * t360 - t35 * (mrSges(7,1) * t159 + mrSges(7,3) * t158) - t105 * (mrSges(6,1) * t159 - mrSges(6,2) * t158) + qJD(6) * t109 - t100 * t101 + qJ(6) * t29 - pkin(5) * t31 + (t158 * t26 + t159 * t27) * mrSges(7,2) + t375; t159 * t101 - t174 * t109 + 0.2e1 * (t4 / 0.2e1 + t35 * t358 + t27 * t357) * m(7) + t31;];
tauc  = t1(:);

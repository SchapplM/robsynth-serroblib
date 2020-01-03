% Calculate vector of inverse dynamics joint torques for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:03
% EndTime: 2019-12-31 19:43:32
% DurationCPUTime: 16.60s
% Computational Cost: add. (3519->583), mult. (8057->766), div. (0->0), fcn. (5236->8), ass. (0->281)
t399 = Ifges(5,1) + Ifges(4,1);
t193 = cos(qJ(2));
t179 = t193 * qJDD(1);
t190 = sin(qJ(2));
t274 = qJD(1) * qJD(2);
t263 = t190 * t274;
t148 = -t179 + t263;
t331 = t148 / 0.2e1;
t149 = qJDD(1) * t190 + t193 * t274;
t187 = sin(pkin(8));
t188 = cos(pkin(8));
t112 = qJDD(2) * t187 + t149 * t188;
t333 = t112 / 0.2e1;
t397 = Ifges(4,5) + Ifges(5,4);
t404 = t331 * t397 + t399 * t333;
t284 = qJD(1) * t190;
t267 = t187 * t284;
t140 = -t188 * qJD(2) + t267;
t403 = -t140 / 0.2e1;
t283 = qJD(1) * t193;
t166 = qJD(5) + t283;
t330 = -t166 / 0.2e1;
t266 = t188 * t284;
t141 = qJD(2) * t187 + t266;
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t224 = t140 * t189 + t141 * t192;
t338 = -t224 / 0.2e1;
t74 = t140 * t192 - t141 * t189;
t340 = -t74 / 0.2e1;
t401 = Ifges(6,5) * t338 + Ifges(6,6) * t340 + Ifges(6,3) * t330;
t111 = -t188 * qJDD(2) + t149 * t187;
t334 = t111 / 0.2e1;
t400 = -mrSges(5,2) - mrSges(4,3);
t398 = -Ifges(4,4) + Ifges(5,5);
t380 = Ifges(4,6) - Ifges(5,6);
t379 = -Ifges(4,3) - Ifges(5,2);
t309 = Ifges(5,5) * t187;
t311 = Ifges(4,4) * t187;
t396 = t188 * t399 + t309 - t311;
t395 = m(6) * pkin(7) + mrSges(6,3) + t400;
t240 = t188 * mrSges(5,1) + t187 * mrSges(5,3);
t242 = mrSges(4,1) * t188 - mrSges(4,2) * t187;
t144 = t187 * t192 - t188 * t189;
t223 = t187 * t189 + t188 * t192;
t389 = t223 * mrSges(6,1) + mrSges(6,2) * t144;
t394 = t240 + t242 + t389;
t386 = -t148 / 0.2e1;
t393 = m(4) * qJD(3);
t173 = Ifges(3,4) * t283;
t392 = t187 * (t141 * Ifges(5,5) - Ifges(5,6) * t283 + t140 * Ifges(5,3)) + Ifges(3,1) * t284 + Ifges(3,5) * qJD(2) + t173 + t188 * (t140 * t398 + t141 * t399 - t283 * t397);
t310 = Ifges(4,4) * t188;
t234 = -Ifges(4,2) * t187 + t310;
t239 = mrSges(5,1) * t187 - mrSges(5,3) * t188;
t241 = mrSges(4,1) * t187 + mrSges(4,2) * t188;
t174 = pkin(6) * t284;
t247 = qJD(2) * pkin(2) - qJD(3) - t174;
t291 = t188 * t193;
t293 = t187 * t193;
t213 = qJ(4) * t141 + t247;
t50 = pkin(3) * t140 - t213;
t180 = t190 * qJ(3);
t183 = t193 * pkin(2);
t269 = -pkin(1) - t183;
t222 = t269 - t180;
t134 = t222 * qJD(1);
t175 = pkin(6) * t283;
t156 = qJD(2) * qJ(3) + t175;
t78 = t134 * t188 - t187 * t156;
t54 = pkin(3) * t283 + qJD(4) - t78;
t79 = t187 * t134 + t188 * t156;
t58 = -qJ(4) * t283 + t79;
t391 = t79 * (-mrSges(4,2) * t190 - mrSges(4,3) * t293) + t78 * (mrSges(4,1) * t190 - mrSges(4,3) * t291) + t58 * (-mrSges(5,2) * t293 + mrSges(5,3) * t190) + t54 * (-mrSges(5,1) * t190 + mrSges(5,2) * t291) - (-t239 * t50 + t241 * t247) * t193 - (t141 * Ifges(4,4) - t140 * Ifges(4,2) - Ifges(4,6) * t283) * t293 / 0.2e1 + (Ifges(4,6) * t190 + t193 * t234) * t403;
t172 = pkin(6) * t179;
t114 = qJDD(2) * qJ(3) + t172 + (qJD(3) - t174) * qJD(2);
t278 = qJD(3) * t190;
t300 = qJDD(1) * pkin(1);
t64 = pkin(2) * t148 - qJ(3) * t149 - qJD(1) * t278 - t300;
t30 = -t187 * t114 + t188 * t64;
t245 = qJDD(4) - t30;
t25 = -pkin(3) * t148 + t245;
t390 = mrSges(5,2) * t25 + t334 * t398 + t404;
t313 = Ifges(3,4) * t190;
t373 = t193 * Ifges(3,2);
t235 = t313 + t373;
t32 = pkin(4) * t283 - pkin(7) * t141 + t54;
t38 = pkin(7) * t140 + t58;
t8 = -t189 * t38 + t192 * t32;
t9 = t189 * t32 + t192 * t38;
t388 = -t8 * mrSges(6,1) + t9 * mrSges(6,2) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t235 / 0.2e1 + t379 * t283 / 0.2e1 + t397 * t141 / 0.2e1 + t380 * t403 + t401;
t19 = qJD(5) * t74 + t111 * t189 + t112 * t192;
t344 = t19 / 0.2e1;
t20 = -qJD(5) * t224 + t111 * t192 - t112 * t189;
t343 = t20 / 0.2e1;
t387 = -m(6) - m(5);
t142 = qJDD(5) - t148;
t332 = t142 / 0.2e1;
t139 = t149 * pkin(6);
t384 = qJD(2) / 0.2e1;
t383 = mrSges(2,2) - mrSges(3,3);
t382 = mrSges(4,2) - mrSges(5,3);
t29 = -mrSges(6,1) * t74 + mrSges(6,2) * t224;
t83 = mrSges(5,1) * t140 - mrSges(5,3) * t141;
t378 = t29 - t83;
t268 = -pkin(6) * t187 - pkin(3);
t204 = -pkin(7) * t291 + (-pkin(4) + t268) * t190;
t229 = pkin(2) * t190 - qJ(3) * t193;
t147 = t229 * qJD(1);
t295 = t147 * t188;
t47 = qJD(1) * t204 - t295;
t129 = t187 * t147;
t170 = qJ(4) * t284;
t292 = t188 * t190;
t217 = -pkin(6) * t292 + pkin(7) * t293;
t55 = qJD(1) * t217 + t129 + t170;
t315 = pkin(7) - qJ(3);
t153 = t315 * t187;
t154 = t315 * t188;
t88 = -t153 * t192 + t154 * t189;
t376 = qJD(3) * t223 + qJD(5) * t88 - t189 * t47 - t192 * t55;
t89 = -t153 * t189 - t154 * t192;
t375 = qJD(3) * t144 - qJD(5) * t89 + t189 * t55 - t192 * t47;
t244 = mrSges(3,1) * t193 - mrSges(3,2) * t190;
t372 = -mrSges(2,1) - t244;
t370 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t140 + mrSges(4,2) * t141 + mrSges(3,3) * t284;
t216 = t144 * t193;
t109 = qJD(1) * t216;
t124 = t144 * qJD(5);
t369 = t109 + t124;
t215 = t223 * t193;
t110 = qJD(1) * t215;
t123 = t223 * qJD(5);
t368 = t110 + t123;
t276 = qJD(4) * t193;
t282 = qJD(2) * t190;
t367 = qJ(4) * t282 - t276;
t138 = -pkin(6) * t263 + t172;
t365 = t138 * t193 + t139 * t190;
t364 = t400 * t190;
t31 = t188 * t114 + t187 * t64;
t363 = -t187 * t30 + t188 * t31;
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t362 = g(1) * t194 + g(2) * t191;
t360 = m(4) - t387;
t256 = qJ(4) * t187 + pkin(2);
t150 = -pkin(3) * t188 - t256;
t324 = pkin(4) * t188;
t358 = t395 * t193 + (-m(6) * (t150 - t324) - m(5) * t150 + m(4) * pkin(2) + t394) * t190;
t308 = Ifges(5,5) * t188;
t230 = Ifges(5,3) * t187 + t308;
t357 = t140 * (Ifges(5,6) * t190 + t193 * t230) + (t190 * t397 + t193 * t396) * t141;
t353 = qJ(3) * t360;
t352 = m(6) * pkin(4) + mrSges(4,1) + mrSges(5,1);
t336 = pkin(3) + pkin(4);
t10 = -pkin(7) * t112 - t148 * t336 + t245;
t21 = t148 * qJ(4) - qJD(1) * t276 + t31;
t11 = pkin(7) * t111 + t21;
t1 = qJD(5) * t8 + t10 * t189 + t11 * t192;
t2 = -qJD(5) * t9 + t10 * t192 - t11 * t189;
t351 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t349 = Ifges(5,5) * t333 + Ifges(5,6) * t331 - t112 * Ifges(4,4) / 0.2e1 + Ifges(4,6) * t386 - t21 * mrSges(5,2) + (Ifges(5,3) + Ifges(4,2)) * t334;
t346 = Ifges(6,4) * t344 + Ifges(6,2) * t343 + Ifges(6,6) * t332;
t345 = Ifges(6,1) * t344 + Ifges(6,4) * t343 + Ifges(6,5) * t332;
t325 = Ifges(6,4) * t224;
t23 = Ifges(6,2) * t74 + Ifges(6,6) * t166 + t325;
t342 = t23 / 0.2e1;
t71 = Ifges(6,4) * t74;
t24 = Ifges(6,1) * t224 + Ifges(6,5) * t166 + t71;
t341 = t24 / 0.2e1;
t339 = t74 / 0.2e1;
t337 = t224 / 0.2e1;
t335 = -t111 / 0.2e1;
t329 = t166 / 0.2e1;
t323 = pkin(6) * t190;
t322 = pkin(7) * t190;
t312 = Ifges(3,4) * t193;
t301 = qJ(4) * t188;
t119 = -qJDD(2) * pkin(2) + qJDD(3) + t139;
t299 = t119 * t190;
t122 = qJD(2) * t229 - t278;
t298 = t122 * t188;
t294 = t187 * t190;
t290 = t190 * t194;
t289 = t191 * t193;
t288 = t193 * t194;
t287 = t194 * t187;
t286 = t183 + t180;
t152 = -pkin(1) - t286;
t104 = pkin(6) * t291 + t187 * t152;
t285 = t194 * pkin(1) + t191 * pkin(6);
t281 = qJD(2) * t193;
t280 = qJD(3) * t187;
t279 = qJD(3) * t188;
t277 = qJD(4) * t187;
t273 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t142;
t272 = pkin(6) * t282;
t265 = qJD(4) * t292;
t260 = -t283 / 0.2e1;
t41 = t111 * mrSges(5,1) - t112 * mrSges(5,3);
t255 = -t274 / 0.2e1;
t254 = t274 / 0.2e1;
t42 = t111 * mrSges(4,1) + t112 * mrSges(4,2);
t63 = -t148 * mrSges(5,1) + t112 * mrSges(5,2);
t160 = pkin(6) * t293;
t103 = t152 * t188 - t160;
t253 = pkin(3) * t291 + qJ(4) * t293 + t286;
t252 = pkin(2) * t288 + qJ(3) * t290 + t285;
t249 = t268 * t190;
t7 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t91 = -qJ(4) * t193 + t104;
t243 = mrSges(3,1) * t190 + mrSges(3,2) * t193;
t117 = t144 * t190;
t118 = t223 * t190;
t238 = mrSges(6,1) * t117 - mrSges(6,2) * t118;
t233 = Ifges(5,4) * t188 + Ifges(5,6) * t187;
t232 = Ifges(3,5) * t193 - Ifges(3,6) * t190;
t231 = Ifges(4,5) * t188 - Ifges(4,6) * t187;
t228 = pkin(3) * t187 - t301;
t45 = -mrSges(6,2) * t166 + mrSges(6,3) * t74;
t46 = mrSges(6,1) * t166 - mrSges(6,3) * t224;
t227 = -t189 * t46 + t192 * t45;
t182 = t193 * pkin(3);
t59 = pkin(4) * t193 + t160 + t182 + (-t152 - t322) * t188;
t77 = pkin(7) * t294 + t91;
t27 = -t189 * t77 + t192 * t59;
t28 = t189 * t59 + t192 * t77;
t125 = t187 * t289 + t188 * t194;
t126 = t188 * t289 - t287;
t226 = t125 * t192 - t126 * t189;
t225 = -t125 * t189 - t126 * t192;
t94 = -pkin(6) * t266 + t129;
t113 = t187 * t122;
t87 = -t188 * t272 + t113;
t221 = pkin(6) + t228;
t220 = pkin(1) * t243;
t219 = -t187 * t336 + t301;
t218 = t190 * (Ifges(3,1) * t193 - t313);
t214 = -pkin(6) + t219;
t127 = -t191 * t188 + t193 * t287;
t205 = -g(1) * t127 - g(2) * t125 - g(3) * t294;
t199 = t193 * (Ifges(5,2) * t190 + t193 * t233);
t198 = t193 * (Ifges(4,3) * t190 + t193 * t231);
t26 = t111 * pkin(3) - t112 * qJ(4) - t141 * qJD(4) + t119;
t184 = t194 * pkin(6);
t157 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t283;
t130 = t188 * t336 + t256;
t128 = t191 * t187 + t188 * t288;
t115 = t221 * t190;
t108 = mrSges(5,1) * t283 + mrSges(5,2) * t141;
t107 = -mrSges(4,1) * t283 - mrSges(4,3) * t141;
t106 = mrSges(4,2) * t283 - mrSges(4,3) * t140;
t105 = -mrSges(5,2) * t140 - mrSges(5,3) * t283;
t97 = t228 * t283 + t175;
t93 = pkin(6) * t267 + t295;
t92 = -t103 + t182;
t90 = t214 * t190;
t86 = t187 * t272 + t298;
t85 = qJD(1) * t249 - t295;
t82 = t170 + t94;
t81 = t221 * t281 - t265;
t80 = t219 * t283 - t175;
t72 = qJD(2) * t249 - t298;
t62 = mrSges(4,1) * t148 - mrSges(4,3) * t112;
t61 = -mrSges(4,2) * t148 - mrSges(4,3) * t111;
t60 = -mrSges(5,2) * t111 + mrSges(5,3) * t148;
t57 = t127 * t189 + t128 * t192;
t56 = t127 * t192 - t128 * t189;
t53 = t214 * t281 + t265;
t51 = t87 + t367;
t49 = qJD(2) * t215 + qJD(5) * t117;
t48 = qJD(2) * t216 - t123 * t190;
t40 = qJD(2) * t217 + t113 + t367;
t39 = qJD(2) * t204 - t298;
t37 = -t140 * t336 + t213;
t14 = pkin(4) * t111 + t26;
t13 = -mrSges(6,2) * t142 + mrSges(6,3) * t20;
t12 = mrSges(6,1) * t142 - mrSges(6,3) * t19;
t6 = -qJD(5) * t28 - t189 * t40 + t192 * t39;
t5 = qJD(5) * t27 + t189 * t39 + t192 * t40;
t3 = [m(4) * (t103 * t30 + t104 * t31 + t78 * t86 + t79 * t87 + (-t247 * t281 + t299) * pkin(6)) + t244 * t300 + t392 * t281 / 0.2e1 + (-mrSges(4,3) * t30 + t390) * t292 + (t232 * t384 + t391) * qJD(2) + (-m(3) * t285 - m(4) * t252 - t57 * mrSges(6,1) - t56 * mrSges(6,2) + t387 * (t128 * pkin(3) + t127 * qJ(4) + t252) + t372 * t194 + t383 * t191 - t352 * t128 + t382 * t127 + t395 * t290) * g(2) + (t26 * t239 - t254 * t373 + Ifges(3,1) * t149 + Ifges(3,4) * t386 + Ifges(3,5) * qJDD(2) + t230 * t334 + t234 * t335 + t396 * t333 + (t231 + t233) * t331) * t190 + (Ifges(6,1) * t49 + Ifges(6,4) * t48) * t337 + (Ifges(6,1) * t118 + Ifges(6,4) * t117) * t344 + (-Ifges(5,6) * t334 - Ifges(4,6) * t335 + pkin(6) * (-qJDD(2) * mrSges(3,2) - mrSges(3,3) * t148) + Ifges(6,3) * t332 + Ifges(6,6) * t343 + Ifges(6,5) * t344 + t312 * t254 - t30 * mrSges(4,1) + t25 * mrSges(5,1) + t31 * mrSges(4,2) - t21 * mrSges(5,3) - t397 * t333 + t379 * t331 + t351) * t193 - (-t380 * t111 + t112 * t397 - t379 * t148) * t193 / 0.2e1 + (-Ifges(6,5) * t337 - Ifges(6,6) * t339 - Ifges(6,3) * t329 + t388) * t282 + t115 * t41 + t87 * t106 + t86 * t107 + t72 * t108 + t103 * t62 + t104 * t61 + t51 * t105 + t81 * t83 + t90 * t7 + t91 * t60 + t92 * t63 + t37 * (-mrSges(6,1) * t48 + mrSges(6,2) * t49) + t53 * t29 + t5 * t45 + t6 * t46 + t27 * t12 + t28 * t13 + (-qJDD(2) * mrSges(3,1) + t42) * t323 + (-t31 * mrSges(4,3) + t349) * t294 + (t1 * t117 - t118 * t2 + t48 * t9 - t49 * t8) * mrSges(6,3) + (t199 + t198) * t255 + m(5) * (t115 * t26 + t21 * t91 + t25 * t92 + t50 * t81 + t51 * t58 + t54 * t72) + m(6) * (t1 * t28 - t14 * t90 + t2 * t27 + t37 * t53 + t5 * t9 + t6 * t8) + qJDD(2) * Ifges(3,6) * t193 + t149 * t312 / 0.2e1 + t370 * pkin(6) * t281 + t218 * t254 + (Ifges(6,5) * t49 + Ifges(6,6) * t48) * t329 + (Ifges(6,5) * t118 + Ifges(6,6) * t117) * t332 + t14 * t238 + t357 * t384 + t235 * t386 + (Ifges(3,4) * t149 - Ifges(3,2) * t148 + t273) * t193 / 0.2e1 + Ifges(2,3) * qJDD(1) + t241 * t299 + (t149 * t323 + t365) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t365 * pkin(6)) + (-t225 * mrSges(6,1) + t226 * mrSges(6,2) + t387 * (-t126 * pkin(3) - qJ(4) * t125 + t184) + t383 * t194 + (-m(3) - m(4)) * t184 + t352 * t126 - t382 * t125 + (m(3) * pkin(1) - m(6) * t269 - (m(6) * t315 + mrSges(6,3)) * t190 + (-m(4) - m(5)) * t222 - t364 - t372) * t191) * g(1) - pkin(1) * (mrSges(3,1) * t148 + mrSges(3,2) * t149) + t49 * t341 + t48 * t342 + t118 * t345 + t117 * t346 + (Ifges(6,4) * t49 + Ifges(6,2) * t48) * t339 + (Ifges(6,4) * t118 + Ifges(6,2) * t117) * t343 - t157 * t272 - t220 * t274; (Ifges(6,5) * t144 - Ifges(6,6) * t223) * t332 + (Ifges(6,4) * t144 - Ifges(6,2) * t223) * t343 + (Ifges(6,1) * t144 - Ifges(6,4) * t223) * t344 + (g(3) * t190 - t1 * t223 - t144 * t2 + t368 * t8 - t369 * t9) * mrSges(6,3) - t223 * t346 + (-pkin(2) * t119 + t363 * qJ(3) + t175 * t247 - t78 * t93 - t79 * t94) * m(4) + (t191 * t358 - t289 * t353) * g(2) + (t194 * t358 - t288 * t353) * g(1) + ((t61 + t60) * qJ(3) + t79 * t393 - Ifges(5,3) * t334 + Ifges(4,2) * t335 + m(5) * (qJ(3) * t21 + qJD(3) * t58) + t380 * t331 - t349) * t188 + (-Ifges(3,2) * t260 - t388 - t401) * t284 + ((t63 - t62) * qJ(3) - t78 * t393 + m(5) * (qJ(3) * t25 + qJD(3) * t54 - qJD(4) * t50) + t390 + t404) * t187 + (t173 + t392) * t260 + (-t357 / 0.2e1 + (t220 - t218 / 0.2e1 + t199 / 0.2e1 + t198 / 0.2e1) * qJD(1) - t391) * qJD(1) - t14 * t389 + (-m(6) * (t253 - t322) - t244 - m(4) * t286 - m(5) * t253 + (-m(6) * t324 - t394) * t193 + t364) * g(3) + (-Ifges(6,5) * t123 - Ifges(6,6) * t124) * t329 + (-Ifges(6,1) * t123 - Ifges(6,4) * t124) * t337 + (-Ifges(6,4) * t123 - Ifges(6,2) * t124) * t339 - Ifges(3,6) * t148 + t130 * t7 - t138 * mrSges(3,2) - t139 * mrSges(3,1) - t109 * t23 / 0.2e1 - t110 * t24 / 0.2e1 - t97 * t83 - t80 * t29 + t88 * t12 + t89 * t13 - pkin(2) * t42 + (Ifges(6,1) * t110 + Ifges(6,4) * t109) * t338 + t309 * t334 + t311 * t335 + (-t82 + t279) * t105 + (-t94 + t279) * t106 + (-t308 + t310) * t333 + t157 * t174 + (t150 * t26 - t50 * t97 - t54 * t85 - t58 * t82) * m(5) + (Ifges(6,4) * t110 + Ifges(6,2) * t109) * t340 + (-t85 + t280) * t108 + t232 * t255 - t26 * t240 - t119 * t242 + (Ifges(6,5) * t110 + Ifges(6,6) * t109) * t330 + t362 * t243 + t363 * mrSges(4,3) + (mrSges(6,1) * t369 - mrSges(6,2) * t368) * t37 - t370 * t175 + t375 * t46 + t376 * t45 + (t1 * t89 - t130 * t14 + t2 * t88 + t376 * t9 + t375 * t8 + (-t80 + t277) * t37) * m(6) + t378 * t277 + Ifges(3,5) * t149 + t150 * t41 - t123 * t341 - t124 * t342 + t144 * t345 + (-t93 - t280) * t107 + Ifges(3,3) * qJDD(2); t74 * t45 - t224 * t46 + (t107 - t108) * t141 + (t105 + t106) * t140 - t7 + t42 + t41 + (-t224 * t8 + t74 * t9 + t14) * m(6) + (t140 * t58 - t141 * t54 + t26) * m(5) + (t140 * t79 + t141 * t78 + t119) * m(4) + (t193 * g(3) - t190 * t362) * t360; t192 * t12 + t189 * t13 - t378 * t141 + t227 * qJD(5) + (t105 + t227) * t283 + t63 + (t1 * t189 - t141 * t37 + t192 * t2 + t205 - t166 * (t189 * t8 - t192 * t9)) * m(6) + (t141 * t50 + t283 * t58 + t205 + t25) * m(5); -t37 * (mrSges(6,1) * t224 + mrSges(6,2) * t74) + (Ifges(6,1) * t74 - t325) * t338 + t23 * t337 + (Ifges(6,5) * t74 - Ifges(6,6) * t224) * t330 - t8 * t45 + t9 * t46 - g(1) * (mrSges(6,1) * t56 - mrSges(6,2) * t57) - g(2) * (mrSges(6,1) * t226 + mrSges(6,2) * t225) - g(3) * t238 + (t224 * t9 + t74 * t8) * mrSges(6,3) + t273 + (-Ifges(6,2) * t224 + t24 + t71) * t340 + t351;];
tau = t3;

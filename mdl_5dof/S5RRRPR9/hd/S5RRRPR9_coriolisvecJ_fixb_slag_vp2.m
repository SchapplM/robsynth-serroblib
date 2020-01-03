% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:15
% EndTime: 2019-12-31 21:22:40
% DurationCPUTime: 11.11s
% Computational Cost: add. (7467->582), mult. (18985->833), div. (0->0), fcn. (13158->8), ass. (0->277)
t239 = sin(pkin(9));
t240 = cos(pkin(9));
t242 = sin(qJ(3));
t245 = cos(qJ(3));
t200 = t239 * t245 + t240 * t242;
t246 = cos(qJ(2));
t250 = t200 * t246;
t165 = qJD(1) * t250;
t186 = t200 * qJD(3);
t293 = -t165 + t186;
t254 = t239 * t242 - t240 * t245;
t251 = t246 * t254;
t166 = qJD(1) * t251;
t187 = t254 * qJD(3);
t292 = -t166 + t187;
t315 = -qJ(4) - pkin(7);
t268 = qJD(3) * t315;
t281 = qJD(4) * t245;
t183 = t242 * t268 + t281;
t184 = -qJD(4) * t242 + t245 * t268;
t121 = -t183 * t239 + t240 * t184;
t243 = sin(qJ(2));
t264 = pkin(2) * t243 - pkin(7) * t246;
t208 = t264 * qJD(1);
t289 = qJD(1) * t243;
t275 = t242 * t289;
t159 = pkin(6) * t275 + t245 * t208;
t294 = t245 * t246;
t252 = pkin(3) * t243 - qJ(4) * t294;
t127 = qJD(1) * t252 + t159;
t190 = t242 * t208;
t295 = t243 * t245;
t296 = t242 * t246;
t144 = t190 + (-pkin(6) * t295 - qJ(4) * t296) * qJD(1);
t70 = t240 * t127 - t144 * t239;
t391 = t121 - t70;
t122 = t240 * t183 + t239 * t184;
t71 = t239 * t127 + t240 * t144;
t390 = t122 - t71;
t241 = sin(qJ(5));
t244 = cos(qJ(5));
t288 = qJD(1) * t246;
t229 = qJD(3) - t288;
t286 = qJD(2) * t245;
t206 = -t275 + t286;
t274 = t245 * t289;
t207 = qJD(2) * t242 + t274;
t143 = t206 * t239 + t207 * t240;
t379 = pkin(8) * t143;
t211 = -pkin(2) * t246 - t243 * pkin(7) - pkin(1);
t194 = t211 * qJD(1);
t236 = pkin(6) * t288;
t218 = qJD(2) * pkin(7) + t236;
t148 = t245 * t194 - t218 * t242;
t113 = -qJ(4) * t207 + t148;
t102 = pkin(3) * t229 + t113;
t149 = t194 * t242 + t218 * t245;
t114 = qJ(4) * t206 + t149;
t107 = t239 * t114;
t50 = t240 * t102 - t107;
t38 = pkin(4) * t229 - t379 + t50;
t267 = t240 * t206 - t207 * t239;
t370 = pkin(8) * t267;
t298 = t240 * t114;
t51 = t239 * t102 + t298;
t39 = t51 + t370;
t10 = t241 * t38 + t244 * t39;
t287 = qJD(2) * t243;
t269 = qJD(1) * t287;
t280 = qJD(2) * qJD(3);
t283 = qJD(3) * t243;
t285 = qJD(2) * t246;
t157 = t245 * t280 + (-t242 * t283 + t245 * t285) * qJD(1);
t282 = qJD(3) * t245;
t374 = t242 * t285 + t243 * t282;
t158 = -qJD(1) * t374 - t242 * t280;
t100 = -t157 * t239 + t158 * t240;
t101 = t157 * t240 + t158 * t239;
t372 = -t143 * t241 + t244 * t267;
t28 = qJD(5) * t372 + t100 * t241 + t101 * t244;
t80 = t143 * t244 + t241 * t267;
t29 = -qJD(5) * t80 + t100 * t244 - t101 * t241;
t279 = Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * t269;
t327 = Ifges(6,4) * t80;
t219 = qJD(5) + t229;
t334 = -t219 / 0.2e1;
t73 = Ifges(6,4) * t372;
t34 = Ifges(6,1) * t80 + Ifges(6,5) * t219 + t73;
t354 = -t80 / 0.2e1;
t356 = -t372 / 0.2e1;
t209 = t264 * qJD(2);
t195 = qJD(1) * t209;
t266 = pkin(6) * t269;
t89 = -qJD(3) * t149 + t245 * t195 + t242 * t266;
t49 = pkin(3) * t269 - qJ(4) * t157 - qJD(4) * t207 + t89;
t284 = qJD(3) * t242;
t88 = t194 * t282 + t242 * t195 - t218 * t284 - t245 * t266;
t54 = qJ(4) * t158 + qJD(4) * t206 + t88;
t15 = -t239 * t54 + t240 * t49;
t11 = pkin(4) * t269 - pkin(8) * t101 + t15;
t16 = t239 * t49 + t240 * t54;
t12 = pkin(8) * t100 + t16;
t9 = -t241 * t39 + t244 * t38;
t2 = qJD(5) * t9 + t11 * t241 + t12 * t244;
t3 = -qJD(5) * t10 + t11 * t244 - t12 * t241;
t369 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t217 = -qJD(2) * pkin(2) + pkin(6) * t289;
t156 = -pkin(3) * t206 + qJD(4) + t217;
t90 = -pkin(4) * t267 + t156;
t389 = t279 + t369 + (Ifges(6,5) * t372 - Ifges(6,6) * t80) * t334 + (t10 * t80 + t372 * t9) * mrSges(6,3) + (-Ifges(6,2) * t80 + t34 + t73) * t356 - t90 * (mrSges(6,1) * t80 + mrSges(6,2) * t372) + (Ifges(6,1) * t372 - t327) * t354;
t388 = -pkin(4) * t289 + t292 * pkin(8) + t391;
t387 = t293 * pkin(8) - t390;
t33 = Ifges(6,2) * t372 + Ifges(6,6) * t219 + t327;
t385 = t33 / 0.2e1;
t271 = Ifges(3,5) * qJD(2) / 0.2e1;
t383 = Ifges(4,3) + Ifges(5,3);
t215 = t315 * t242;
t216 = t315 * t245;
t152 = t240 * t215 + t216 * t239;
t125 = -pkin(8) * t200 + t152;
t153 = t239 * t215 - t240 * t216;
t126 = -pkin(8) * t254 + t153;
t61 = t125 * t241 + t126 * t244;
t378 = -qJD(5) * t61 + t387 * t241 + t388 * t244;
t60 = t125 * t244 - t126 * t241;
t377 = qJD(5) * t60 + t388 * t241 - t387 * t244;
t376 = Ifges(5,4) * t143;
t324 = pkin(3) * t242;
t197 = t288 * t324 + t236;
t375 = pkin(3) * t284 + t293 * pkin(4) - t197;
t234 = Ifges(3,4) * t288;
t306 = t207 * Ifges(4,4);
t129 = t206 * Ifges(4,2) + t229 * Ifges(4,6) + t306;
t198 = Ifges(4,4) * t206;
t130 = t207 * Ifges(4,1) + t229 * Ifges(4,5) + t198;
t255 = t148 * t245 + t149 * t242;
t312 = Ifges(4,4) * t245;
t259 = -Ifges(4,2) * t242 + t312;
t313 = Ifges(4,4) * t242;
t261 = Ifges(4,1) * t245 - t313;
t262 = mrSges(4,1) * t242 + mrSges(4,2) * t245;
t310 = Ifges(4,6) * t242;
t311 = Ifges(4,5) * t245;
t329 = t245 / 0.2e1;
t330 = -t242 / 0.2e1;
t331 = t229 / 0.2e1;
t335 = t207 / 0.2e1;
t247 = -t255 * mrSges(4,3) + t130 * t329 + t206 * t259 / 0.2e1 + t261 * t335 + t217 * t262 + (-t310 + t311) * t331 + t129 * t330;
t373 = t247 + Ifges(3,1) * t289 / 0.2e1 + t234 / 0.2e1 + t271;
t68 = Ifges(5,2) * t267 + Ifges(5,6) * t229 + t376;
t371 = t68 / 0.2e1;
t345 = -t267 / 0.2e1;
t270 = -Ifges(3,6) * qJD(2) / 0.2e1;
t368 = Ifges(5,4) * t267;
t232 = pkin(3) * t240 + pkin(4);
t325 = pkin(3) * t239;
t182 = t232 * t241 + t244 * t325;
t56 = -t113 * t239 - t298;
t42 = t56 - t370;
t57 = t240 * t113 - t107;
t43 = t57 - t379;
t367 = -t182 * qJD(5) + t241 * t43 - t244 * t42;
t181 = t232 * t244 - t241 * t325;
t366 = t181 * qJD(5) - t241 * t42 - t244 * t43;
t363 = Ifges(4,5) * t157 + Ifges(4,6) * t158;
t231 = pkin(6) * t294;
t168 = t242 * t211 + t231;
t362 = Ifges(5,5) * t101 + Ifges(5,6) * t100 + t269 * t383 + t363;
t361 = qJD(1) * pkin(1) * mrSges(3,2);
t277 = Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t314 = Ifges(3,4) * t243;
t360 = -t277 * t229 - t148 * mrSges(4,1) - t50 * mrSges(5,1) - t9 * mrSges(6,1) - t207 * Ifges(4,5) - t206 * Ifges(4,6) - t270 + (Ifges(3,2) * t246 + t314) * qJD(1) / 0.2e1 - t219 * Ifges(6,3) - t80 * Ifges(6,5) - t372 * Ifges(6,6) - t143 * Ifges(5,5) - t267 * Ifges(5,6) + t10 * mrSges(6,2) + t149 * mrSges(4,2) + t51 * mrSges(5,2) - t383 * t331;
t359 = -t89 * mrSges(4,1) - t15 * mrSges(5,1) + t88 * mrSges(4,2) + t16 * mrSges(5,2);
t358 = t28 / 0.2e1;
t357 = t29 / 0.2e1;
t355 = t372 / 0.2e1;
t353 = t80 / 0.2e1;
t352 = pkin(1) * mrSges(3,1);
t349 = t100 / 0.2e1;
t348 = t101 / 0.2e1;
t176 = t200 * t243;
t177 = t254 * t243;
t117 = -t176 * t244 + t177 * t241;
t347 = t117 / 0.2e1;
t118 = -t176 * t241 - t177 * t244;
t346 = t118 / 0.2e1;
t344 = t267 / 0.2e1;
t343 = -t143 / 0.2e1;
t342 = t143 / 0.2e1;
t341 = t157 / 0.2e1;
t340 = t158 / 0.2e1;
t339 = -t176 / 0.2e1;
t338 = -t177 / 0.2e1;
t337 = -t206 / 0.2e1;
t336 = -t207 / 0.2e1;
t333 = t219 / 0.2e1;
t332 = -t229 / 0.2e1;
t326 = pkin(3) * t207;
t323 = pkin(6) * t242;
t290 = t245 * t209 + t287 * t323;
t72 = -t243 * t281 + t252 * qJD(2) + (-t231 + (qJ(4) * t243 - t211) * t242) * qJD(3) + t290;
t291 = t242 * t209 + t211 * t282;
t81 = (-qJD(2) * pkin(6) - qJ(4) * qJD(3)) * t295 + (-qJD(4) * t243 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t246) * t242 + t291;
t36 = t239 * t72 + t240 * t81;
t105 = -t165 * t244 + t166 * t241;
t136 = t200 * t244 - t241 * t254;
t75 = -qJD(5) * t136 - t186 * t244 + t187 * t241;
t303 = t105 - t75;
t106 = -t165 * t241 - t166 * t244;
t135 = -t200 * t241 - t244 * t254;
t74 = qJD(5) * t135 - t186 * t241 - t187 * t244;
t302 = t106 - t74;
t299 = qJD(2) * mrSges(3,2);
t297 = t242 * t243;
t202 = t245 * t211;
t145 = -qJ(4) * t295 + t202 + (-pkin(3) - t323) * t246;
t151 = -qJ(4) * t297 + t168;
t84 = t239 * t145 + t240 * t151;
t210 = pkin(3) * t297 + t243 * pkin(6);
t164 = pkin(3) * t374 + pkin(6) * t285;
t233 = -pkin(3) * t245 - pkin(2);
t8 = -t29 * mrSges(6,1) + t28 * mrSges(6,2);
t46 = -t100 * mrSges(5,1) + t101 * mrSges(5,2);
t139 = -pkin(3) * t158 + qJD(2) * t236;
t35 = -t239 * t81 + t240 * t72;
t83 = t240 * t145 - t239 * t151;
t265 = m(4) * t217 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t206 + mrSges(4,2) * t207 + mrSges(3,3) * t289;
t263 = mrSges(4,1) * t245 - mrSges(4,2) * t242;
t260 = Ifges(4,1) * t242 + t312;
t258 = Ifges(4,2) * t245 + t313;
t257 = Ifges(4,5) * t242 + Ifges(4,6) * t245;
t58 = -pkin(4) * t246 + t177 * pkin(8) + t83;
t62 = -pkin(8) * t176 + t84;
t30 = -t241 * t62 + t244 * t58;
t31 = t241 * t58 + t244 * t62;
t256 = -t242 * t89 + t245 * t88;
t213 = mrSges(3,3) * t288 - t299;
t169 = pkin(4) * t254 + t233;
t167 = -pkin(6) * t296 + t202;
t163 = mrSges(4,1) * t229 - mrSges(4,3) * t207;
t162 = -mrSges(4,2) * t229 + mrSges(4,3) * t206;
t160 = -pkin(6) * t274 + t190;
t146 = pkin(4) * t176 + t210;
t138 = -mrSges(4,2) * t269 + mrSges(4,3) * t158;
t137 = mrSges(4,1) * t269 - mrSges(4,3) * t157;
t120 = -qJD(2) * t251 - t186 * t243;
t119 = -qJD(2) * t250 + t254 * t283;
t116 = mrSges(5,1) * t229 - mrSges(5,3) * t143;
t115 = -mrSges(5,2) * t229 + mrSges(5,3) * t267;
t112 = -qJD(3) * t168 + t290;
t111 = (-t243 * t286 - t246 * t284) * pkin(6) + t291;
t104 = pkin(4) * t143 + t326;
t103 = -mrSges(4,1) * t158 + mrSges(4,2) * t157;
t92 = t157 * Ifges(4,1) + t158 * Ifges(4,4) + Ifges(4,5) * t269;
t91 = t157 * Ifges(4,4) + t158 * Ifges(4,2) + Ifges(4,6) * t269;
t87 = mrSges(5,1) * t269 - mrSges(5,3) * t101;
t86 = -mrSges(5,2) * t269 + mrSges(5,3) * t100;
t85 = -pkin(4) * t119 + t164;
t82 = -mrSges(5,1) * t267 + mrSges(5,2) * t143;
t69 = Ifges(5,1) * t143 + Ifges(5,5) * t229 + t368;
t64 = mrSges(6,1) * t219 - mrSges(6,3) * t80;
t63 = -mrSges(6,2) * t219 + mrSges(6,3) * t372;
t59 = -pkin(4) * t100 + t139;
t45 = t101 * Ifges(5,1) + t100 * Ifges(5,4) + Ifges(5,5) * t269;
t44 = t101 * Ifges(5,4) + t100 * Ifges(5,2) + Ifges(5,6) * t269;
t41 = -qJD(5) * t118 + t119 * t244 - t120 * t241;
t40 = qJD(5) * t117 + t119 * t241 + t120 * t244;
t37 = -mrSges(6,1) * t372 + mrSges(6,2) * t80;
t24 = pkin(8) * t119 + t36;
t23 = pkin(4) * t287 - pkin(8) * t120 + t35;
t22 = -mrSges(6,2) * t269 + mrSges(6,3) * t29;
t21 = mrSges(6,1) * t269 - mrSges(6,3) * t28;
t7 = t28 * Ifges(6,1) + t29 * Ifges(6,4) + Ifges(6,5) * t269;
t6 = t28 * Ifges(6,4) + t29 * Ifges(6,2) + Ifges(6,6) * t269;
t5 = -qJD(5) * t31 + t23 * t244 - t24 * t241;
t4 = qJD(5) * t30 + t23 * t241 + t24 * t244;
t1 = [-(t279 + t362 + t363) * t246 / 0.2e1 + (t92 * t329 + pkin(6) * t103 + t91 * t330 + t261 * t341 + t259 * t340 + (-t242 * t88 - t245 * t89) * mrSges(4,3) + (t258 * t337 + t260 * t336 + t217 * t263 + t257 * t332 - t245 * t129 / 0.2e1 + t130 * t330 + (t148 * t242 - t149 * t245) * mrSges(4,3)) * qJD(3) + (((t311 / 0.2e1 - t310 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,4)) * t243 + Ifges(5,5) * t338 + Ifges(5,6) * t339 + Ifges(6,5) * t346 + Ifges(6,6) * t347 - 0.2e1 * t352) * qJD(1) - pkin(6) * t213 + t270 - t360) * qJD(2)) * t243 + m(4) * (t149 * t111 + t148 * t112 + t89 * t167 + t88 * t168) + (-Ifges(5,1) * t177 - Ifges(5,4) * t176) * t348 + (-Ifges(5,4) * t177 - Ifges(5,2) * t176) * t349 + (t119 * t51 - t120 * t50 + t15 * t177 - t16 * t176) * mrSges(5,3) + t139 * (mrSges(5,1) * t176 - mrSges(5,2) * t177) + t120 * t69 / 0.2e1 + t36 * t115 + t35 * t116 + t59 * (-mrSges(6,1) * t117 + mrSges(6,2) * t118) + t84 * t86 + t83 * t87 + t90 * (-mrSges(6,1) * t41 + mrSges(6,2) * t40) + t85 * t37 + t4 * t63 + t5 * t64 + t40 * t34 / 0.2e1 + t30 * t21 + t31 * t22 + (t10 * t41 + t117 * t2 - t118 * t3 - t40 * t9) * mrSges(6,3) + t41 * t385 + (Ifges(6,1) * t118 + Ifges(6,4) * t117) * t358 + (-Ifges(5,5) * t348 - Ifges(6,5) * t358 - Ifges(5,6) * t349 - Ifges(6,6) * t357 + (0.3e1 / 0.2e1 * Ifges(3,4) * t285 + (-Ifges(6,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(6) + t262) * pkin(6) - t277) * t287) * qJD(1) + t359 - t369) * t246 + (Ifges(6,4) * t118 + Ifges(6,2) * t117) * t357 + (t265 * pkin(6) + t271 - 0.2e1 * t361 + t373) * t285 + t146 * t8 + t156 * (-mrSges(5,1) * t119 + mrSges(5,2) * t120) + t111 * t162 + t112 * t163 + t164 * t82 + t167 * t137 + t168 * t138 + t210 * t46 + (Ifges(5,5) * t120 + Ifges(5,6) * t119) * t331 + (Ifges(6,5) * t40 + Ifges(6,6) * t41) * t333 + t45 * t338 + t44 * t339 + (Ifges(5,1) * t120 + Ifges(5,4) * t119) * t342 + (Ifges(5,4) * t120 + Ifges(5,2) * t119) * t344 + t7 * t346 + t6 * t347 + (Ifges(6,1) * t40 + Ifges(6,4) * t41) * t353 + (Ifges(6,4) * t40 + Ifges(6,2) * t41) * t355 + t119 * t371 + m(6) * (t10 * t4 + t146 * t59 + t2 * t31 + t3 * t30 + t5 * t9 + t85 * t90) + m(5) * (t139 * t210 + t15 * t83 + t156 * t164 + t16 * t84 + t35 * t50 + t36 * t51); m(5) * (t121 * t50 + t122 * t51 + t139 * t233 + t15 * t152 + t153 * t16) + (-t15 * t200 - t16 * t254 + t292 * t50 - t293 * t51) * mrSges(5,3) + t139 * (mrSges(5,1) * t254 + mrSges(5,2) * t200) + (Ifges(5,1) * t200 - Ifges(5,4) * t254) * t348 + (Ifges(5,4) * t200 - Ifges(5,2) * t254) * t349 - t254 * t44 / 0.2e1 + (t165 / 0.2e1 - t186 / 0.2e1) * t68 + (-Ifges(5,5) * t166 - Ifges(5,6) * t165) * t332 + (-Ifges(5,1) * t166 - Ifges(5,4) * t165) * t343 + (-Ifges(5,4) * t166 - Ifges(5,2) * t165) * t345 + (t166 / 0.2e1 - t187 / 0.2e1) * t69 + t256 * mrSges(4,3) + (mrSges(5,1) * t293 - mrSges(5,2) * t292) * t156 - m(4) * (t148 * t159 + t149 * t160) + (-Ifges(5,5) * t187 - Ifges(5,6) * t186) * t331 + (-Ifges(5,1) * t187 - Ifges(5,4) * t186) * t342 + (-Ifges(5,4) * t187 - Ifges(5,2) * t186) * t344 + t390 * t115 + t391 * t116 - pkin(2) * t103 + t60 * t21 + t61 * t22 - m(5) * (t156 * t197 + t50 * t70 + t51 * t71) + ((m(5) * t156 + t82) * t324 + t247) * qJD(3) + (m(4) * t256 - t242 * t137 + t245 * t138 + (-m(4) * t255 - t242 * t162 - t245 * t163) * qJD(3)) * pkin(7) + t377 * t63 + t378 * t64 + (t10 * t377 + t169 * t59 + t2 * t61 + t3 * t60 + t375 * t90 + t378 * t9) * m(6) + t375 * t37 + ((t361 - t234 / 0.2e1 + t271 + ((-m(4) * pkin(2) - mrSges(3,1) - t263) * qJD(2) - t265) * pkin(6) - t373) * t246 + ((t352 + t314 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t246) * qJD(1) + (t213 + t299) * pkin(6) + t270 + t360) * t243 + (Ifges(5,5) * t200 + Ifges(6,5) * t136 - Ifges(5,6) * t254 + Ifges(6,6) * t135 + t257) * t287 / 0.2e1) * qJD(1) + (-t106 / 0.2e1 + t74 / 0.2e1) * t34 + (mrSges(6,1) * t303 - mrSges(6,2) * t302) * t90 + (-t10 * t303 + t135 * t2 - t136 * t3 + t302 * t9) * mrSges(6,3) + t135 * t6 / 0.2e1 + t136 * t7 / 0.2e1 + t59 * (-mrSges(6,1) * t135 + mrSges(6,2) * t136) + t152 * t87 + t153 * t86 - t160 * t162 - t159 * t163 + t169 * t8 - t197 * t82 + t200 * t45 / 0.2e1 + t233 * t46 + t91 * t329 + (Ifges(6,5) * t74 + Ifges(6,6) * t75) * t333 + (Ifges(6,5) * t106 + Ifges(6,6) * t105) * t334 + t258 * t340 + t260 * t341 + (Ifges(6,1) * t74 + Ifges(6,4) * t75) * t353 + (Ifges(6,1) * t106 + Ifges(6,4) * t105) * t354 + (Ifges(6,4) * t74 + Ifges(6,2) * t75) * t355 + (Ifges(6,4) * t106 + Ifges(6,2) * t105) * t356 + (Ifges(6,4) * t136 + Ifges(6,2) * t135) * t357 + (Ifges(6,1) * t136 + Ifges(6,4) * t135) * t358 + t242 * t92 / 0.2e1 + (-t105 / 0.2e1 + t75 / 0.2e1) * t33; (-t207 * t82 + t239 * t86 + t240 * t87) * pkin(3) + t389 + (-Ifges(4,2) * t207 + t130 + t198) * t337 + (-Ifges(5,2) * t143 + t368 + t69) * t345 + (t143 * t51 + t267 * t50) * mrSges(5,3) + (Ifges(4,5) * t206 + Ifges(5,5) * t267 - Ifges(4,6) * t207 - Ifges(5,6) * t143) * t332 - t156 * (mrSges(5,1) * t143 + mrSges(5,2) * t267) + t143 * t371 - t359 + t366 * t63 + t367 * t64 + (t10 * t366 - t104 * t90 + t181 * t3 + t182 * t2 + t367 * t9) * m(6) + t80 * t385 - t57 * t115 - t56 * t116 - t104 * t37 + ((t15 * t240 + t16 * t239) * pkin(3) - t156 * t326 - t50 * t56 - t51 * t57) * m(5) + (t148 * t206 + t149 * t207) * mrSges(4,3) + (Ifges(5,1) * t267 - t376) * t343 + t362 - t148 * t162 + t149 * t163 + t181 * t21 + t182 * t22 - t217 * (mrSges(4,1) * t207 + mrSges(4,2) * t206) + t129 * t335 + (Ifges(4,1) * t206 - t306) * t336; -t267 * t115 + t143 * t116 - t372 * t63 + t80 * t64 + t46 + t8 + (-t10 * t372 + t80 * t9 + t59) * m(6) + (t143 * t50 - t267 * t51 + t139) * m(5); t10 * t64 + t33 * t353 - t9 * t63 + t389;];
tauc = t1(:);

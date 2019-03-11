% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:11
% EndTime: 2019-03-09 09:34:34
% DurationCPUTime: 11.49s
% Computational Cost: add. (10123->625), mult. (23436->867), div. (0->0), fcn. (15695->8), ass. (0->281)
t266 = sin(qJ(2));
t255 = t266 * qJD(1);
t268 = cos(qJ(5));
t262 = cos(pkin(10));
t303 = t262 * t255;
t261 = sin(pkin(10));
t265 = sin(qJ(5));
t326 = t261 * t265;
t172 = -t255 * t326 + t268 * t303;
t316 = qJD(5) * t268;
t317 = qJD(5) * t265;
t196 = -t261 * t317 + t262 * t316;
t322 = t196 + t172;
t283 = t268 * t261 + t265 * t262;
t277 = t283 * t266;
t173 = qJD(1) * t277;
t197 = t283 * qJD(5);
t424 = t197 + t173;
t250 = pkin(2) * t255;
t269 = cos(qJ(2));
t287 = -qJ(3) * t269 + qJ(4) * t266;
t187 = qJD(1) * t287 + t250;
t320 = qJD(1) * t269;
t251 = pkin(7) * t320;
t220 = pkin(3) * t320 + t251;
t134 = -t261 * t187 + t262 * t220;
t325 = t261 * t266;
t281 = pkin(4) * t269 - pkin(8) * t325;
t105 = qJD(1) * t281 + t134;
t135 = t262 * t187 + t261 * t220;
t118 = pkin(8) * t303 + t135;
t263 = -pkin(2) - qJ(4);
t349 = -pkin(8) + t263;
t223 = t349 * t261;
t224 = t349 * t262;
t395 = -qJD(4) * t283 - t265 * t105 - t268 * t118 - t223 * t317 + t224 * t316;
t157 = t268 * t223 + t265 * t224;
t282 = -t262 * t268 + t326;
t394 = qJD(4) * t282 - qJD(5) * t157 - t268 * t105 + t265 * t118;
t264 = sin(qJ(6));
t267 = cos(qJ(6));
t315 = t262 * qJD(2);
t211 = t261 * t320 - t315;
t212 = -t261 * qJD(2) - t262 * t320;
t297 = t211 * t265 + t268 * t212;
t300 = -t266 * qJ(3) - pkin(1);
t210 = t263 * t269 + t300;
t174 = t210 * qJD(1);
t249 = pkin(7) * t255;
t219 = -pkin(3) * t255 - t249;
t389 = -t219 + qJD(3);
t180 = qJD(2) * t263 + t389;
t113 = -t174 * t261 + t262 * t180;
t85 = pkin(4) * t255 + pkin(8) * t211 + t113;
t114 = t262 * t174 + t261 * t180;
t92 = pkin(8) * t212 + t114;
t42 = t265 * t85 + t268 * t92;
t36 = pkin(9) * t297 + t42;
t332 = t267 * t36;
t244 = t255 + qJD(5);
t143 = t211 * t268 - t212 * t265;
t41 = -t265 * t92 + t268 * t85;
t35 = pkin(9) * t143 + t41;
t34 = pkin(5) * t244 + t35;
t10 = t264 * t34 + t332;
t148 = t264 * t282 - t267 * t283;
t276 = t281 * qJD(2);
t313 = qJD(1) * qJD(2);
t302 = t266 * t313;
t243 = pkin(2) * t302;
t314 = t266 * qJD(3);
t272 = qJD(2) * t287 - qJD(4) * t269 - t314;
t142 = qJD(1) * t272 + t243;
t380 = pkin(3) + pkin(7);
t305 = t380 * t269;
t186 = (qJD(1) * t305 - qJD(4)) * qJD(2);
t97 = -t261 * t142 + t262 * t186;
t66 = qJD(1) * t276 + t97;
t295 = t262 * t302;
t98 = t262 * t142 + t261 * t186;
t82 = pkin(8) * t295 + t98;
t16 = -qJD(5) * t42 - t265 * t82 + t268 * t66;
t301 = t269 * t313;
t275 = qJD(2) * t277;
t99 = qJD(1) * t275 + qJD(5) * t297;
t11 = pkin(5) * t301 - t99 * pkin(9) + t16;
t319 = qJD(2) * t266;
t274 = t282 * t319;
t100 = -qJD(1) * t274 + qJD(5) * t143;
t15 = t265 * t66 + t268 * t82 + t85 * t316 - t317 * t92;
t12 = pkin(9) * t100 + t15;
t334 = t264 * t36;
t9 = t267 * t34 - t334;
t2 = qJD(6) * t9 + t11 * t264 + t12 * t267;
t3 = -qJD(6) * t10 + t11 * t267 - t12 * t264;
t273 = qJD(6) * t148 - t196 * t264 - t267 * t197;
t298 = t172 * t264 + t267 * t173;
t330 = t298 - t273;
t391 = -t264 * t283 - t267 * t282;
t106 = t172 * t267 - t173 * t264;
t74 = qJD(6) * t391 + t196 * t267 - t197 * t264;
t420 = t106 + t74;
t423 = -t10 * t420 + t148 * t2 - t391 * t3 + t330 * t9;
t422 = -pkin(5) * t320 + pkin(9) * t424 + t394;
t421 = pkin(9) * t322 - t395;
t415 = t143 * t264 + t267 * t297;
t32 = qJD(6) * t415 + t100 * t264 + t267 * t99;
t69 = t143 * t267 - t264 * t297;
t33 = qJD(6) * t69 + t100 * t267 - t264 * t99;
t312 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t301;
t354 = Ifges(7,4) * t69;
t237 = qJD(6) + t244;
t362 = -t237 / 0.2e1;
t383 = t69 / 0.2e1;
t385 = -t415 / 0.2e1;
t67 = Ifges(7,4) * t415;
t39 = -Ifges(7,1) * t69 + Ifges(7,5) * t237 + t67;
t404 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t260 = qJD(2) * qJ(3);
t192 = qJD(4) + t260 + t220;
t153 = -pkin(4) * t212 + t192;
t86 = -pkin(5) * t297 + t153;
t419 = t312 + t404 + (Ifges(7,5) * t415 + Ifges(7,6) * t69) * t362 + (-t10 * t69 + t415 * t9) * mrSges(7,3) + (Ifges(7,2) * t69 + t39 + t67) * t385 - t86 * (-mrSges(7,1) * t69 + mrSges(7,2) * t415) + (Ifges(7,1) * t415 + t354) * t383;
t407 = -qJD(1) / 0.2e1;
t418 = -mrSges(3,1) + mrSges(4,2);
t304 = -pkin(4) * t262 - pkin(3);
t177 = t255 * t304 - t249;
t417 = qJD(3) - t177;
t416 = -t15 * t283 + t16 * t282 - t322 * t42 + t41 * t424;
t38 = Ifges(7,2) * t415 + Ifges(7,6) * t237 - t354;
t412 = t38 / 0.2e1;
t406 = -qJD(2) / 0.2e1;
t405 = qJD(2) / 0.2e1;
t156 = -t223 * t265 + t268 * t224;
t121 = pkin(9) * t282 + t156;
t122 = -pkin(9) * t283 + t157;
t55 = t121 * t264 + t122 * t267;
t403 = -qJD(6) * t55 + t421 * t264 + t267 * t422;
t54 = t121 * t267 - t122 * t264;
t402 = qJD(6) * t54 + t264 * t422 - t421 * t267;
t235 = t380 * t266;
t216 = t262 * t235;
t128 = t266 * pkin(4) + t216 + (pkin(8) * t269 - t210) * t261;
t152 = t262 * t210 + t261 * t235;
t323 = t262 * t269;
t132 = -pkin(8) * t323 + t152;
t58 = t265 * t128 + t268 * t132;
t392 = pkin(5) * t322 + t417;
t390 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t227 = -pkin(2) * t269 + t300;
t206 = t227 * qJD(1);
t231 = -t251 - t260;
t233 = -mrSges(4,1) * t320 - qJD(2) * mrSges(4,3);
t286 = t113 * t261 - t114 * t262;
t346 = Ifges(5,4) * t261;
t290 = Ifges(5,2) * t262 + t346;
t345 = Ifges(5,4) * t262;
t347 = Ifges(5,1) * t261;
t291 = t345 + t347;
t348 = mrSges(5,2) * t261;
t292 = mrSges(5,1) * t262 - t348;
t356 = t262 / 0.2e1;
t357 = t261 / 0.2e1;
t388 = t286 * mrSges(5,3) - (m(4) * t231 + qJD(2) * mrSges(3,2) - mrSges(3,3) * t320 + t233) * pkin(7) + t192 * t292 - t231 * mrSges(4,1) - Ifges(3,6) * t406 - Ifges(4,5) * t405 + t206 * mrSges(4,2) + t211 * t291 / 0.2e1 - t212 * t290 / 0.2e1 - (-t211 * Ifges(5,1) + t212 * Ifges(5,4) + Ifges(5,5) * t255) * t357 - (-t211 * Ifges(5,4) + t212 * Ifges(5,2) + Ifges(5,6) * t255) * t356 + ((-Ifges(3,2) - Ifges(4,3)) * t269 + (-Ifges(3,4) - Ifges(4,6)) * t266) * t407;
t387 = t32 / 0.2e1;
t386 = t33 / 0.2e1;
t384 = t415 / 0.2e1;
t382 = -t69 / 0.2e1;
t381 = t99 / 0.2e1;
t379 = pkin(1) * mrSges(3,1);
t378 = pkin(1) * mrSges(3,2);
t376 = t100 / 0.2e1;
t182 = t282 * t269;
t183 = t283 * t269;
t125 = t182 * t267 + t183 * t264;
t375 = t125 / 0.2e1;
t126 = t182 * t264 - t183 * t267;
t374 = t126 / 0.2e1;
t373 = -t297 / 0.2e1;
t372 = t297 / 0.2e1;
t371 = t143 / 0.2e1;
t370 = -t143 / 0.2e1;
t369 = t148 / 0.2e1;
t368 = t391 / 0.2e1;
t367 = -(Ifges(5,6) * t269 + t266 * t290) * t313 / 0.2e1;
t366 = t182 / 0.2e1;
t365 = -t183 / 0.2e1;
t364 = -t283 / 0.2e1;
t363 = -t282 / 0.2e1;
t361 = t237 / 0.2e1;
t360 = -t244 / 0.2e1;
t359 = t244 / 0.2e1;
t358 = -t261 / 0.2e1;
t257 = t269 * pkin(7);
t344 = Ifges(6,4) * t143;
t343 = Ifges(5,5) * t261;
t342 = Ifges(4,6) * t269;
t341 = Ifges(5,6) * t262;
t333 = t266 * mrSges(4,3);
t245 = t261 * pkin(4) + qJ(3);
t253 = pkin(2) * t319;
t161 = t253 + t272;
t222 = qJD(2) * t305;
t117 = t262 * t161 + t261 * t222;
t236 = t269 * pkin(3) + t257;
t318 = qJD(2) * t269;
t311 = Ifges(6,5) * t99 + Ifges(6,6) * t100 + Ifges(6,3) * t301;
t310 = 0.3e1 / 0.2e1 * Ifges(3,4) + 0.3e1 / 0.2e1 * Ifges(4,6);
t309 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t308 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t155 = -mrSges(5,1) * t212 - mrSges(5,2) * t211;
t80 = -mrSges(6,1) * t297 - mrSges(6,2) * t143;
t307 = -t233 + t155 + t80;
t306 = m(4) * pkin(7) + mrSges(4,1);
t195 = pkin(4) * t323 + t236;
t8 = -t33 * mrSges(7,1) + t32 * mrSges(7,2);
t49 = -t100 * mrSges(6,1) + t99 * mrSges(6,2);
t57 = t268 * t128 - t132 * t265;
t116 = -t261 * t161 + t262 * t222;
t221 = t380 * t319;
t226 = -qJD(2) * pkin(2) + qJD(3) + t249;
t294 = m(4) * t226 + (mrSges(4,1) + mrSges(3,3)) * t255 + t418 * qJD(2);
t288 = t261 * t98 + t262 * t97;
t50 = pkin(5) * t266 + pkin(9) * t183 + t57;
t51 = pkin(9) * t182 + t58;
t23 = -t264 * t51 + t267 * t50;
t24 = t264 * t50 + t267 * t51;
t279 = -t343 / 0.2e1 - t341 / 0.2e1;
t178 = (-pkin(7) + t304) * t319;
t278 = -qJ(3) * t318 - t314;
t102 = pkin(8) * t266 * t315 + t117;
t96 = t116 + t276;
t25 = t268 * t102 + t128 * t316 - t132 * t317 + t265 * t96;
t176 = -mrSges(5,1) * t295 + t302 * t348;
t259 = qJD(2) * qJD(3);
t158 = qJD(1) * t178 + t259;
t26 = -qJD(5) * t58 - t265 * t102 + t268 * t96;
t248 = Ifges(3,4) * t320;
t270 = t113 * mrSges(5,1) + t226 * mrSges(4,1) + t41 * mrSges(6,1) + t9 * mrSges(7,1) + t212 * Ifges(5,6) - t211 * Ifges(5,5) + Ifges(3,5) * t405 + t248 / 0.2e1 + Ifges(4,4) * t406 + (-Ifges(4,2) * t266 - t342) * t407 + t237 * Ifges(7,3) - t69 * Ifges(7,5) + t415 * Ifges(7,6) + t244 * Ifges(6,3) - t143 * Ifges(6,5) + t297 * Ifges(6,6) - t10 * mrSges(7,2) - t114 * mrSges(5,2) - t206 * mrSges(4,3) - t42 * mrSges(6,2) + (Ifges(5,3) + Ifges(3,1)) * t255 / 0.2e1;
t225 = pkin(7) * t302 - t259;
t217 = (mrSges(4,2) * t269 - t333) * qJD(1);
t190 = t253 + t278;
t189 = (mrSges(5,3) * t262 * t266 - mrSges(5,2) * t269) * t313;
t188 = (mrSges(5,1) * t269 - mrSges(5,3) * t325) * t313;
t185 = -qJD(1) * t221 + t259;
t175 = qJD(1) * t278 + t243;
t171 = pkin(5) * t283 + t245;
t170 = mrSges(5,1) * t255 + mrSges(5,3) * t211;
t169 = -mrSges(5,2) * t255 + mrSges(5,3) * t212;
t160 = (Ifges(5,5) * t269 + t266 * t291) * t313;
t151 = -t261 * t210 + t216;
t141 = Ifges(6,4) * t297;
t140 = -pkin(5) * t182 + t195;
t131 = t197 * t269 - t274;
t130 = qJD(5) * t182 + t275;
t124 = mrSges(6,1) * t244 + mrSges(6,3) * t143;
t123 = -mrSges(6,2) * t244 + mrSges(6,3) * t297;
t91 = -pkin(5) * t131 + t178;
t84 = -mrSges(6,2) * t301 + t100 * mrSges(6,3);
t83 = mrSges(6,1) * t301 - t99 * mrSges(6,3);
t63 = -Ifges(6,1) * t143 + Ifges(6,5) * t244 + t141;
t62 = Ifges(6,2) * t297 + Ifges(6,6) * t244 - t344;
t60 = mrSges(7,1) * t237 + mrSges(7,3) * t69;
t59 = -mrSges(7,2) * t237 + mrSges(7,3) * t415;
t56 = -pkin(5) * t100 + t158;
t47 = Ifges(6,1) * t99 + Ifges(6,4) * t100 + Ifges(6,5) * t301;
t46 = Ifges(6,4) * t99 + Ifges(6,2) * t100 + Ifges(6,6) * t301;
t44 = -qJD(6) * t126 - t130 * t264 + t131 * t267;
t43 = qJD(6) * t125 + t130 * t267 + t131 * t264;
t40 = -mrSges(7,1) * t415 - mrSges(7,2) * t69;
t28 = -mrSges(7,2) * t301 + t33 * mrSges(7,3);
t27 = mrSges(7,1) * t301 - t32 * mrSges(7,3);
t20 = pkin(9) * t131 + t25;
t19 = pkin(5) * t318 - t130 * pkin(9) + t26;
t14 = t267 * t35 - t334;
t13 = -t264 * t35 - t332;
t7 = Ifges(7,1) * t32 + Ifges(7,4) * t33 + Ifges(7,5) * t301;
t6 = Ifges(7,4) * t32 + Ifges(7,2) * t33 + Ifges(7,6) * t301;
t5 = -qJD(6) * t24 + t19 * t267 - t20 * t264;
t4 = qJD(6) * t23 + t19 * t264 + t20 * t267;
t1 = [t44 * t412 + (t175 * mrSges(4,2) + t185 * t292 - t225 * mrSges(4,1) + t160 * t358 + t262 * t367 + (t261 * t97 - t262 * t98) * mrSges(5,3) + (t294 * pkin(7) + qJD(2) * t309 + t270) * qJD(2) + (-t227 * mrSges(4,3) - 0.2e1 * t378 + Ifges(7,5) * t374 + Ifges(7,6) * t375 + Ifges(6,5) * t365 + Ifges(6,6) * t366 + (t279 + t310) * t269) * t313) * t269 + (t312 + t311) * t266 / 0.2e1 + (t10 * t44 + t125 * t2 - t126 * t3 - t43 * t9) * mrSges(7,3) + m(5) * (t113 * t116 + t114 * t117 + t151 * t97 + t152 * t98 + t185 * t236 - t192 * t221) + (-t130 * t41 + t131 * t42 + t15 * t182 + t16 * t183) * mrSges(6,3) + (-Ifges(6,4) * t183 + Ifges(6,2) * t182) * t376 + (-Ifges(6,1) * t183 + Ifges(6,4) * t182) * t381 + t158 * (-mrSges(6,1) * t182 - mrSges(6,2) * t183) + (Ifges(7,4) * t126 + Ifges(7,2) * t125) * t386 + (Ifges(7,1) * t126 + Ifges(7,4) * t125) * t387 + t236 * t176 + t190 * t217 - t221 * t155 + t195 * t49 + t151 * t188 + t152 * t189 + t178 * t80 + t117 * t169 + t116 * t170 + t153 * (-mrSges(6,1) * t131 + mrSges(6,2) * t130) + t140 * t8 + t131 * t62 / 0.2e1 + t56 * (-mrSges(7,1) * t125 + mrSges(7,2) * t126) + t130 * t63 / 0.2e1 + t25 * t123 + t26 * t124 + (Ifges(7,4) * t43 + Ifges(7,2) * t44) * t384 + (Ifges(6,4) * t130 + Ifges(6,2) * t131) * t372 + t7 * t374 + t6 * t375 + (Ifges(7,1) * t43 + Ifges(7,4) * t44) * t382 + (Ifges(6,1) * t130 + Ifges(6,4) * t131) * t370 + (Ifges(6,5) * t130 + Ifges(6,6) * t131) * t359 + (Ifges(7,5) * t43 + Ifges(7,6) * t44) * t361 + t47 * t365 + t46 * t366 + t86 * (-mrSges(7,1) * t44 + mrSges(7,2) * t43) + t91 * t40 + t57 * t83 + t58 * t84 + t4 * t59 + t5 * t60 + m(7) * (t10 * t4 + t140 * t56 + t2 * t24 + t23 * t3 + t5 * t9 + t86 * t91) + m(6) * (t15 * t58 + t153 * t178 + t158 * t195 + t16 * t57 + t25 * t42 + t26 * t41) + (-t98 * mrSges(5,2) + t97 * mrSges(5,1) + Ifges(7,6) * t386 + Ifges(7,5) * t387 + Ifges(6,6) * t376 + Ifges(6,5) * t381 + (qJD(2) * t308 - t388) * qJD(2) + (-t227 * mrSges(4,2) - 0.2e1 * t379 + (0.3e1 / 0.2e1 * t343 + 0.3e1 / 0.2e1 * t341 - t310) * t266 + (0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - t262 ^ 2 * Ifges(5,2) / 0.2e1 + Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + t306 * pkin(7) + (-t345 - t347 / 0.2e1) * t261) * t269) * t313 + t390 + t404) * t266 - t175 * t333 + t23 * t27 + t24 * t28 + t43 * t39 / 0.2e1 + m(4) * (t175 * t227 + t206 * t190 - t225 * t257); (Ifges(7,4) * t273 - Ifges(7,2) * t74) * t384 + (Ifges(7,1) * t273 - Ifges(7,4) * t74) * t382 + (Ifges(7,5) * t273 - Ifges(7,6) * t74) * t361 + (-t106 / 0.2e1 - t74 / 0.2e1) * t38 + (((t379 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1 + t279) * t266) * qJD(1) + (pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,1) + (Ifges(5,1) * t262 - t346) * t357 + (-Ifges(5,2) * t261 + t345) * t356 + t308) * qJD(2) + t388) * t266 + ((-pkin(2) * mrSges(4,1) + Ifges(5,5) * t356 + Ifges(6,5) * t363 + Ifges(7,5) * t368 + Ifges(5,6) * t358 + Ifges(6,6) * t364 + Ifges(7,6) * t369 + t309) * qJD(2) + ((-m(4) * pkin(2) + t418) * qJD(2) - t294) * pkin(7) - t248 / 0.2e1 + (-t342 / 0.2e1 + t378) * qJD(1) + (-Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t255 - t270) * t269) * qJD(1) + (mrSges(6,1) * t322 - mrSges(6,2) * t424) * t153 + (-m(4) * t206 - t217) * (-qJ(3) * t320 + t250) + (-t298 / 0.2e1 + t273 / 0.2e1) * t39 + (Ifges(7,1) * t298 + Ifges(7,4) * t106) * t383 + (Ifges(7,4) * t298 + Ifges(7,2) * t106) * t385 + (Ifges(7,5) * t298 + Ifges(7,6) * t106) * t362 + t423 * mrSges(7,3) + t392 * t40 + (mrSges(7,1) * t420 - mrSges(7,2) * t330) * t86 + m(4) * (-t225 * qJ(3) - t231 * qJD(3)) + (t185 * mrSges(5,2) + t160 / 0.2e1 - t97 * mrSges(5,3) - qJD(4) * t170 + t263 * t188) * t262 + t158 * (mrSges(6,1) * t283 - mrSges(6,2) * t282) + (-Ifges(6,4) * t282 - Ifges(6,2) * t283) * t376 + (-Ifges(6,1) * t282 - Ifges(6,4) * t283) * t381 + (Ifges(7,4) * t391 + Ifges(7,2) * t148) * t386 + (Ifges(7,1) * t391 + Ifges(7,4) * t148) * t387 + t56 * (-mrSges(7,1) * t148 + mrSges(7,2) * t391) + (-t113 * t134 - t114 * t135 + qJ(3) * t185 + t288 * t263 + (-t113 * t262 - t114 * t261) * qJD(4) + t389 * t192) * m(5) + (-t196 / 0.2e1 - t172 / 0.2e1) * t62 + (-Ifges(6,1) * t197 - Ifges(6,4) * t196) * t370 + (-Ifges(6,4) * t197 - Ifges(6,2) * t196) * t372 + (-Ifges(6,5) * t197 - Ifges(6,6) * t196) * t359 + (-t197 / 0.2e1 - t173 / 0.2e1) * t63 + (t15 * t157 + t153 * t417 + t156 * t16 + t158 * t245 + t394 * t41 + t395 * t42) * m(6) + t416 * mrSges(6,3) + t245 * t49 - t225 * mrSges(4,3) - t219 * t155 + qJ(3) * t176 - t177 * t80 - t135 * t169 - t134 * t170 + t171 * t8 + t156 * t83 + t157 * t84 + (Ifges(6,4) * t173 + Ifges(6,2) * t172) * t373 + t7 * t368 + t6 * t369 + (Ifges(6,1) * t173 + Ifges(6,4) * t172) * t371 + (Ifges(6,5) * t173 + Ifges(6,6) * t172) * t360 + t47 * t363 + t46 * t364 + (t185 * mrSges(5,1) - t98 * mrSges(5,3) - qJD(4) * t169 + t263 * t189 + t367) * t261 + t402 * t59 + (t10 * t402 + t171 * t56 + t2 * t55 + t3 * t54 + t392 * t86 + t403 * t9) * m(7) + t403 * t60 + t307 * qJD(3) + t54 * t27 + t55 * t28 + t394 * t124 + t395 * t123; t391 * t27 - t148 * t28 + t262 * t188 + t261 * t189 + t283 * t84 - t282 * t83 - t330 * t60 + t420 * t59 - t424 * t124 + t322 * t123 + (t169 * t262 - t170 * t261 + t217) * t255 + (t306 * t320 - t307 - t40) * qJD(2) - m(4) * (-qJD(2) * t231 - t206 * t255) + (-qJD(2) * t86 - t423) * m(7) + (-qJD(2) * t153 - t416) * m(6) + (-qJD(2) * t192 - t255 * t286 + t288) * m(5); -t297 * t123 - t143 * t124 - t212 * t169 - t211 * t170 - t415 * t59 - t69 * t60 + t176 + t49 + t8 + (-t10 * t415 - t69 * t9 + t56) * m(7) + (-t143 * t41 - t297 * t42 + t158) * m(6) + (-t113 * t211 - t114 * t212 + t185) * m(5); t390 - m(7) * (t10 * t14 + t13 * t9) + t311 - t69 * t412 + (-t143 * t42 + t297 * t41) * mrSges(6,3) + (t143 * t40 + t264 * t28 + t267 * t27 + (-t264 * t60 + t267 * t59) * qJD(6) + (t143 * t86 + t2 * t264 + t267 * t3 + (t10 * t267 - t264 * t9) * qJD(6)) * m(7)) * pkin(5) + (Ifges(6,5) * t297 + Ifges(6,6) * t143) * t360 - t153 * (-mrSges(6,1) * t143 + mrSges(6,2) * t297) + (Ifges(6,2) * t143 + t141 + t63) * t373 + (Ifges(6,1) * t297 + t344) * t371 - t41 * t123 + t42 * t124 + t62 * t370 - t13 * t60 - t14 * t59 + t419; t10 * t60 + t38 * t382 - t9 * t59 + t419;];
tauc  = t1(:);

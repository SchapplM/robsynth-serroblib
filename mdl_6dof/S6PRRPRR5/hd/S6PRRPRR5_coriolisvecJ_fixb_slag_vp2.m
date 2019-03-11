% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:24
% EndTime: 2019-03-08 22:16:48
% DurationCPUTime: 12.81s
% Computational Cost: add. (10001->633), mult. (25376->912), div. (0->0), fcn. (19252->12), ass. (0->286)
t262 = sin(qJ(3));
t266 = cos(qJ(3));
t280 = pkin(3) * t262 - qJ(4) * t266;
t201 = qJD(3) * t280 - qJD(4) * t262;
t256 = sin(pkin(12));
t258 = cos(pkin(12));
t263 = sin(qJ(2));
t309 = qJD(3) * t262;
t301 = pkin(8) * t309;
t257 = sin(pkin(6));
t314 = qJD(1) * t257;
t267 = cos(qJ(2));
t319 = t266 * t267;
t387 = t258 * t201 + t256 * t301 - (-t256 * t319 + t258 * t263) * t314;
t418 = t256 * t201 - (t256 * t263 + t258 * t319) * t314;
t321 = t258 * t266;
t277 = pkin(4) * t262 - pkin(9) * t321;
t273 = t277 * qJD(3);
t417 = t273 + t387;
t323 = t258 * t262;
t326 = t256 * t266;
t416 = (-pkin(8) * t323 - pkin(9) * t326) * qJD(3) + t418;
t261 = sin(qJ(5));
t265 = cos(qJ(5));
t224 = t256 * t265 + t258 * t261;
t275 = t224 * t266;
t191 = qJD(2) * t275;
t204 = t224 * qJD(5);
t318 = -t191 + t204;
t322 = t258 * t265;
t278 = t256 * t261 - t322;
t274 = t278 * t266;
t192 = qJD(2) * t274;
t203 = t278 * qJD(5);
t317 = -t192 + t203;
t299 = t263 * t314;
t229 = qJD(2) * pkin(8) + t299;
t214 = t262 * t229;
t259 = cos(pkin(6));
t320 = t259 * t266;
t297 = qJD(1) * t320;
t184 = -t214 + t297;
t226 = t280 * qJD(2);
t131 = -t184 * t256 + t258 * t226;
t108 = qJD(2) * t277 + t131;
t132 = t258 * t184 + t256 * t226;
t310 = qJD(2) * t266;
t296 = t256 * t310;
t115 = -pkin(9) * t296 + t132;
t351 = pkin(9) + qJ(4);
t232 = t351 * t256;
t233 = t351 * t258;
t306 = qJD(5) * t265;
t405 = qJD(4) * t322 - t265 * t115 - t232 * t306 + (-qJD(4) * t256 - qJD(5) * t233 - t108) * t261;
t176 = -t261 * t232 + t265 * t233;
t404 = -t224 * qJD(4) - qJD(5) * t176 - t265 * t108 + t115 * t261;
t311 = qJD(2) * t262;
t415 = -pkin(5) * t311 + pkin(10) * t317 + t404;
t414 = pkin(10) * t318 - t405;
t231 = -pkin(3) * t266 - qJ(4) * t262 - pkin(2);
t217 = t258 * t231;
t159 = -pkin(9) * t323 + t217 + (-pkin(8) * t256 - pkin(4)) * t266;
t188 = pkin(8) * t321 + t256 * t231;
t327 = t256 * t262;
t172 = -pkin(9) * t327 + t188;
t95 = t261 * t159 + t265 * t172;
t390 = -qJD(5) * t95 - t416 * t261 + t265 * t417;
t307 = qJD(5) * t261;
t389 = t159 * t306 - t172 * t307 + t261 * t417 + t416 * t265;
t413 = -t310 / 0.2e1;
t271 = qJD(3) * t274;
t136 = -t204 * t262 - t271;
t412 = pkin(5) * t309 - pkin(10) * t136 + t390;
t272 = qJD(3) * t275;
t137 = t203 * t262 - t272;
t411 = -pkin(10) * t137 - t389;
t264 = cos(qJ(6));
t250 = qJD(5) - t310;
t218 = t258 * qJD(3) - t256 * t311;
t305 = t256 * qJD(3);
t219 = t258 * t311 + t305;
t154 = t218 * t261 + t219 * t265;
t313 = qJD(1) * t262;
t185 = t266 * t229 + t259 * t313;
t181 = qJD(3) * qJ(4) + t185;
t298 = t267 * t314;
t186 = qJD(2) * t231 - t298;
t109 = -t181 * t256 + t258 * t186;
t80 = -pkin(4) * t310 - pkin(9) * t219 + t109;
t110 = t258 * t181 + t256 * t186;
t91 = pkin(9) * t218 + t110;
t44 = -t261 * t91 + t265 * t80;
t33 = -pkin(10) * t154 + t44;
t32 = pkin(5) * t250 + t33;
t260 = sin(qJ(6));
t289 = t265 * t218 - t219 * t261;
t45 = t261 * t80 + t265 * t91;
t34 = pkin(10) * t289 + t45;
t337 = t260 * t34;
t13 = t264 * t32 - t337;
t336 = t264 * t34;
t14 = t260 * t32 + t336;
t105 = -qJD(2) * t271 + qJD(5) * t289;
t106 = -qJD(2) * t272 - qJD(5) * t154;
t402 = -t154 * t260 + t264 * t289;
t29 = qJD(6) * t402 + t105 * t264 + t106 * t260;
t304 = qJD(2) * qJD(3);
t292 = t262 * t304;
t87 = t154 * t264 + t260 * t289;
t30 = -qJD(6) * t87 - t105 * t260 + t106 * t264;
t302 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t292;
t358 = Ifges(7,4) * t87;
t243 = qJD(6) + t250;
t365 = -t243 / 0.2e1;
t378 = -t87 / 0.2e1;
t380 = -t402 / 0.2e1;
t81 = Ifges(7,4) * t402;
t39 = Ifges(7,1) * t87 + Ifges(7,5) * t243 + t81;
t312 = qJD(2) * t257;
t294 = t267 * t312;
t286 = t266 * t294;
t316 = qJD(1) * t286 + qJD(3) * t297;
t139 = (qJD(4) - t214) * qJD(3) + t316;
t164 = (t201 + t299) * qJD(2);
t77 = -t139 * t256 + t258 * t164;
t63 = qJD(2) * t273 + t77;
t291 = t266 * t304;
t285 = t256 * t291;
t78 = t258 * t139 + t256 * t164;
t69 = -pkin(9) * t285 + t78;
t18 = -qJD(5) * t45 - t261 * t69 + t265 * t63;
t10 = pkin(5) * t292 - pkin(10) * t105 + t18;
t17 = t261 * t63 + t265 * t69 + t80 * t306 - t307 * t91;
t12 = pkin(10) * t106 + t17;
t2 = qJD(6) * t13 + t10 * t260 + t12 * t264;
t3 = -qJD(6) * t14 + t10 * t264 - t12 * t260;
t395 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t177 = -qJD(3) * pkin(3) + qJD(4) - t184;
t133 = -pkin(4) * t218 + t177;
t79 = -pkin(5) * t289 + t133;
t410 = t302 + t395 + (Ifges(7,5) * t402 - Ifges(7,6) * t87) * t365 + (t13 * t402 + t14 * t87) * mrSges(7,3) + (-Ifges(7,2) * t87 + t39 + t81) * t380 - t79 * (mrSges(7,1) * t87 + mrSges(7,2) * t402) + (Ifges(7,1) * t402 - t358) * t378;
t409 = Ifges(4,5) / 0.2e1;
t408 = Ifges(4,4) * t413;
t175 = -t265 * t232 - t233 * t261;
t134 = -pkin(10) * t224 + t175;
t135 = -pkin(10) * t278 + t176;
t66 = t134 * t264 - t135 * t260;
t407 = qJD(6) * t66 + t260 * t415 - t414 * t264;
t67 = t134 * t260 + t135 * t264;
t406 = -qJD(6) * t67 + t414 * t260 + t264 * t415;
t162 = pkin(4) * t296 + t185;
t403 = pkin(5) * t318 - t162;
t38 = Ifges(7,2) * t402 + Ifges(7,6) * t243 + t358;
t400 = t38 / 0.2e1;
t399 = t311 / 0.2e1;
t293 = -Ifges(4,6) * qJD(3) / 0.2e1;
t398 = qJD(3) * t409;
t198 = t278 * t262;
t94 = t265 * t159 - t172 * t261;
t64 = -pkin(5) * t266 + pkin(10) * t198 + t94;
t197 = t224 * t262;
t68 = -pkin(10) * t197 + t95;
t36 = t260 * t64 + t264 * t68;
t397 = -qJD(6) * t36 + t260 * t411 + t264 * t412;
t35 = -t260 * t68 + t264 * t64;
t396 = qJD(6) * t35 + t260 * t412 - t264 * t411;
t386 = -t258 * t301 + t418;
t385 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t11 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t350 = mrSges(5,2) * t258;
t194 = mrSges(5,1) * t285 + t291 * t350;
t55 = -t106 * mrSges(6,1) + t105 * mrSges(6,2);
t384 = t11 + t194 + t55;
t230 = -qJD(2) * pkin(2) - t298;
t345 = Ifges(5,2) * t256;
t348 = Ifges(5,4) * t258;
t281 = -t345 + t348;
t349 = Ifges(5,4) * t256;
t282 = Ifges(5,1) * t258 - t349;
t283 = mrSges(5,1) * t256 + t350;
t360 = t258 / 0.2e1;
t361 = -t256 / 0.2e1;
t383 = (t109 * t258 + t110 * t256) * mrSges(5,3) - t177 * t283 - t230 * mrSges(4,2) - Ifges(4,1) * t399 + t408 - t398 + t184 * mrSges(4,3) - t218 * t281 / 0.2e1 - t219 * t282 / 0.2e1 - (t219 * Ifges(5,4) + t218 * Ifges(5,2) - Ifges(5,6) * t310) * t361 - (t219 * Ifges(5,1) + t218 * Ifges(5,4) - Ifges(5,5) * t310) * t360;
t382 = t29 / 0.2e1;
t381 = t30 / 0.2e1;
t379 = t402 / 0.2e1;
t377 = t87 / 0.2e1;
t375 = t105 / 0.2e1;
t374 = t106 / 0.2e1;
t126 = -t197 * t264 + t198 * t260;
t373 = t126 / 0.2e1;
t127 = -t197 * t260 - t198 * t264;
t372 = t127 / 0.2e1;
t371 = -t289 / 0.2e1;
t370 = t289 / 0.2e1;
t369 = -t154 / 0.2e1;
t368 = t154 / 0.2e1;
t367 = -t197 / 0.2e1;
t366 = -t198 / 0.2e1;
t364 = t243 / 0.2e1;
t363 = -t250 / 0.2e1;
t362 = t250 / 0.2e1;
t255 = t262 * pkin(8);
t347 = Ifges(6,4) * t154;
t346 = Ifges(5,5) * t258;
t344 = Ifges(5,6) * t256;
t116 = -t191 * t264 + t192 * t260;
t156 = t224 * t264 - t260 * t278;
t89 = -qJD(6) * t156 + t203 * t260 - t204 * t264;
t335 = t116 - t89;
t117 = -t191 * t260 - t192 * t264;
t155 = -t224 * t260 - t264 * t278;
t88 = qJD(6) * t155 - t203 * t264 - t204 * t260;
t334 = t117 - t88;
t142 = -t229 * t309 + t316;
t331 = t142 * t266;
t308 = qJD(3) * t266;
t143 = t229 * t308 + (qJD(3) * t259 + t294) * t313;
t325 = t257 * t263;
t205 = t262 * t325 - t320;
t330 = t143 * t205;
t324 = t257 * t267;
t315 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t218 - mrSges(5,2) * t219 - mrSges(4,3) * t311;
t212 = t266 * pkin(4) * t305 + pkin(8) * t308;
t227 = pkin(4) * t327 + t255;
t303 = t143 * t255;
t300 = Ifges(6,5) * t105 + Ifges(6,6) * t106 + Ifges(6,3) * t292;
t252 = -pkin(4) * t258 - pkin(3);
t295 = t263 * t312;
t46 = -mrSges(7,1) * t402 + mrSges(7,2) * t87;
t90 = -mrSges(6,1) * t289 + mrSges(6,2) * t154;
t288 = t90 + t46 - t315;
t287 = t262 * t298;
t284 = -mrSges(4,1) * t266 + mrSges(4,2) * t262;
t206 = t259 * t262 + t266 * t325;
t168 = -t206 * t256 - t258 * t324;
t169 = t206 * t258 - t256 * t324;
t96 = t168 * t265 - t169 * t261;
t97 = t168 * t261 + t169 * t265;
t48 = -t260 * t97 + t264 * t96;
t49 = t260 * t96 + t264 * t97;
t118 = pkin(4) * t285 + t143;
t269 = t109 * mrSges(5,1) + t13 * mrSges(7,1) + t230 * mrSges(4,1) + t44 * mrSges(6,1) + Ifges(5,3) * t413 + t219 * Ifges(5,5) + t218 * Ifges(5,6) + t293 - (Ifges(4,4) * t262 + t266 * Ifges(4,2)) * qJD(2) / 0.2e1 + t243 * Ifges(7,3) + t87 * Ifges(7,5) + t402 * Ifges(7,6) + t250 * Ifges(6,3) + t154 * Ifges(6,5) + t289 * Ifges(6,6) - t110 * mrSges(5,2) - t14 * mrSges(7,2) - t185 * mrSges(4,3) - t45 * mrSges(6,2);
t268 = qJD(2) ^ 2;
t238 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t310;
t225 = t284 * qJD(2);
t211 = (mrSges(4,1) * t262 + mrSges(4,2) * t266) * t304;
t200 = (mrSges(5,1) * t262 - mrSges(5,3) * t321) * t304;
t199 = (-mrSges(5,2) * t262 - mrSges(5,3) * t326) * t304;
t193 = pkin(5) * t278 + t252;
t190 = -mrSges(5,1) * t310 - mrSges(5,3) * t219;
t189 = mrSges(5,2) * t310 + mrSges(5,3) * t218;
t187 = -pkin(8) * t326 + t217;
t174 = (Ifges(5,5) * t262 + t266 * t282) * t304;
t173 = (Ifges(5,6) * t262 + t266 * t281) * t304;
t171 = -qJD(3) * t205 + t286;
t170 = qJD(3) * t206 + t262 * t294;
t161 = pkin(5) * t197 + t227;
t150 = Ifges(6,4) * t289;
t130 = t171 * t258 + t256 * t295;
t129 = -t171 * t256 + t258 * t295;
t124 = mrSges(6,1) * t250 - mrSges(6,3) * t154;
t123 = -mrSges(6,2) * t250 + mrSges(6,3) * t289;
t107 = -pkin(5) * t137 + t212;
t93 = -mrSges(6,2) * t292 + mrSges(6,3) * t106;
t92 = mrSges(6,1) * t292 - mrSges(6,3) * t105;
t75 = t154 * Ifges(6,1) + t250 * Ifges(6,5) + t150;
t74 = Ifges(6,2) * t289 + t250 * Ifges(6,6) + t347;
t71 = mrSges(7,1) * t243 - mrSges(7,3) * t87;
t70 = -mrSges(7,2) * t243 + mrSges(7,3) * t402;
t60 = -pkin(5) * t106 + t118;
t54 = t105 * Ifges(6,1) + t106 * Ifges(6,4) + Ifges(6,5) * t292;
t53 = t105 * Ifges(6,4) + t106 * Ifges(6,2) + Ifges(6,6) * t292;
t52 = -qJD(6) * t127 - t136 * t260 + t137 * t264;
t51 = qJD(6) * t126 + t136 * t264 + t137 * t260;
t41 = -qJD(5) * t97 + t129 * t265 - t130 * t261;
t40 = qJD(5) * t96 + t129 * t261 + t130 * t265;
t24 = -mrSges(7,2) * t292 + mrSges(7,3) * t30;
t23 = mrSges(7,1) * t292 - mrSges(7,3) * t29;
t16 = t264 * t33 - t337;
t15 = -t260 * t33 - t336;
t9 = t29 * Ifges(7,1) + t30 * Ifges(7,4) + Ifges(7,5) * t292;
t8 = t29 * Ifges(7,4) + t30 * Ifges(7,2) + Ifges(7,6) * t292;
t7 = -qJD(6) * t49 - t260 * t40 + t264 * t41;
t6 = qJD(6) * t48 + t260 * t41 + t264 * t40;
t1 = [t171 * t238 + t169 * t199 + t168 * t200 + t130 * t189 + t129 * t190 - t206 * mrSges(4,3) * t292 + t40 * t123 + t41 * t124 + t97 * t93 + t96 * t92 + t6 * t70 + t7 * t71 + t48 * t23 + t49 * t24 + ((-mrSges(3,2) * t268 - t211) * t267 + (-mrSges(3,1) * t268 + qJD(2) * t225) * t263) * t257 + (mrSges(4,3) * t291 + t384) * t205 + t288 * t170 + m(7) * (t13 * t7 + t14 * t6 + t170 * t79 + t2 * t49 + t205 * t60 + t3 * t48) + m(6) * (t118 * t205 + t133 * t170 + t17 * t97 + t18 * t96 + t40 * t45 + t41 * t44) + m(5) * (t109 * t129 + t110 * t130 + t168 * t77 + t169 * t78 + t170 * t177 + t330) + m(4) * (t142 * t206 + t330 - t170 * t184 + t171 * t185 + (t230 - t298) * t295); -(t302 + t300) * t266 / 0.2e1 + (-Ifges(6,4) * t198 - Ifges(6,2) * t197) * t374 + (-t136 * t44 + t137 * t45 - t17 * t197 + t18 * t198) * mrSges(6,3) + (-Ifges(6,1) * t198 - Ifges(6,4) * t197) * t375 + t118 * (mrSges(6,1) * t197 - mrSges(6,2) * t198) + (-m(4) * (t230 * t263 + (-t184 * t262 + t185 * t266) * t267) - t225 * t263 - t238 * t319) * t314 + t386 * t189 + t387 * t190 + (t109 * t387 + t110 * t386 - t177 * t287 + t187 * t77 + t188 * t78 + t303) * m(5) + t389 * t123 + (t118 * t227 + t17 * t95 + t18 * t94 + t389 * t45 + t390 * t44 + (t212 - t287) * t133) * m(6) + t390 * t124 + t396 * t70 + t397 * t71 + (t161 * t60 + t2 * t36 + t3 * t35 + (t107 - t287) * t79 + t396 * t14 + t397 * t13) * m(7) + (-t77 * mrSges(5,1) + t78 * mrSges(5,2) - Ifges(6,5) * t375 - Ifges(7,5) * t382 - Ifges(6,6) * t374 - Ifges(7,6) * t381 + (t398 + (-m(4) * t184 + m(5) * t177 - t315) * pkin(8) + ((0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * t346 + 0.3e1 / 0.2e1 * t344) * t266 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + Ifges(5,1) * t258 ^ 2 / 0.2e1 - Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + (-t348 + t345 / 0.2e1) * t256) * t262) * qJD(2) - t383) * qJD(3) + t385 - t395) * t266 + (t126 * t2 - t127 * t3 - t13 * t51 + t14 * t52) * mrSges(7,3) + ((-m(4) * pkin(2) + t284) * t299 + (Ifges(6,5) * t366 + Ifges(6,6) * t367 + Ifges(7,5) * t372 + Ifges(7,6) * t373 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t346 / 0.2e1 - t344 / 0.2e1) * t262) * t309) * qJD(2) + (pkin(8) * t194 + t173 * t361 + t174 * t360 + (mrSges(4,3) + t283) * t143 + (-t256 * t78 - t258 * t77) * mrSges(5,3) - t288 * t298) * t262 + t227 * t55 - pkin(2) * t211 + t212 * t90 + t188 * t199 + t187 * t200 + t161 * t11 + m(4) * (pkin(8) * t331 + t303) + t136 * t75 / 0.2e1 + t133 * (-mrSges(6,1) * t137 + mrSges(6,2) * t136) + t137 * t74 / 0.2e1 + t60 * (-mrSges(7,1) * t126 + mrSges(7,2) * t127) + t107 * t46 + t94 * t92 + t95 * t93 + t79 * (-mrSges(7,1) * t52 + mrSges(7,2) * t51) + t51 * t39 / 0.2e1 + t35 * t23 + t36 * t24 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t381 + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t382 + mrSges(4,3) * t331 + (t269 + (-m(4) * t185 - t238) * pkin(8) + t293) * t309 + t52 * t400 + (Ifges(6,5) * t136 + Ifges(6,6) * t137) * t362 + (Ifges(7,5) * t51 + Ifges(7,6) * t52) * t364 + t54 * t366 + t53 * t367 + (Ifges(6,1) * t136 + Ifges(6,4) * t137) * t368 + (Ifges(6,4) * t136 + Ifges(6,2) * t137) * t370 + t9 * t372 + t8 * t373 + (Ifges(7,1) * t51 + Ifges(7,4) * t52) * t377 + (Ifges(7,4) * t51 + Ifges(7,2) * t52) * t379; ((t408 + (-t344 + t346) * t310 / 0.2e1 + (t409 + (Ifges(5,1) * t256 + t348) * t360 + (Ifges(5,2) * t258 + t349) * t361) * qJD(3) + t383) * t266 + (-t269 + t293 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t310 + Ifges(4,4) * t399) * t262 + (Ifges(5,5) * t256 + Ifges(6,5) * t224 + Ifges(7,5) * t156 + Ifges(5,6) * t258 - Ifges(6,6) * t278 + Ifges(7,6) * t155) * t309 / 0.2e1) * qJD(2) + (t174 / 0.2e1 + t143 * mrSges(5,2) - t77 * mrSges(5,3) - qJD(4) * t190 - qJ(4) * t200) * t256 + (-t143 * mrSges(5,1) + t173 / 0.2e1 + t78 * mrSges(5,3) + qJD(4) * t189 + qJ(4) * t199) * t258 + (-Ifges(6,5) * t203 - Ifges(6,6) * t204) * t362 + (-Ifges(6,1) * t203 - Ifges(6,4) * t204) * t368 + (-Ifges(6,4) * t203 - Ifges(6,2) * t204) * t370 + (-t116 / 0.2e1 + t89 / 0.2e1) * t38 + t403 * t46 + (-t204 / 0.2e1 + t191 / 0.2e1) * t74 + (-t203 / 0.2e1 + t192 / 0.2e1) * t75 + (-Ifges(6,5) * t192 - Ifges(6,6) * t191) * t363 + (-Ifges(6,1) * t192 - Ifges(6,4) * t191) * t369 + (-Ifges(6,4) * t192 - Ifges(6,2) * t191) * t371 + t406 * t71 + t407 * t70 + (t13 * t406 + t14 * t407 + t193 * t60 + t2 * t67 + t3 * t66 + t403 * t79) * m(7) + (-t109 * t131 - t110 * t132 - t177 * t185 - pkin(3) * t143 + (-t109 * t256 + t110 * t258) * qJD(4) + (-t256 * t77 + t258 * t78) * qJ(4)) * m(5) - t278 * t53 / 0.2e1 + t118 * (mrSges(6,1) * t278 + mrSges(6,2) * t224) + (-t17 * t278 - t18 * t224 + t317 * t44 - t318 * t45) * mrSges(6,3) + (Ifges(6,4) * t224 - Ifges(6,2) * t278) * t374 + (Ifges(6,1) * t224 - Ifges(6,4) * t278) * t375 + (-t117 / 0.2e1 + t88 / 0.2e1) * t39 + t404 * t124 + (t118 * t252 - t133 * t162 + t17 * t176 + t175 * t18 + t404 * t44 + t405 * t45) * m(6) + t405 * t123 + t252 * t55 - t184 * t238 + t224 * t54 / 0.2e1 + t193 * t11 - pkin(3) * t194 - t132 * t189 - t131 * t190 + (t13 * t334 - t14 * t335 + t155 * t2 - t156 * t3) * mrSges(7,3) + t175 * t92 + t176 * t93 - t162 * t90 + (mrSges(7,1) * t335 - mrSges(7,2) * t334) * t79 + t155 * t8 / 0.2e1 + t60 * (-mrSges(7,1) * t155 + mrSges(7,2) * t156) + t156 * t9 / 0.2e1 - t142 * mrSges(4,2) - t143 * mrSges(4,1) + (mrSges(6,1) * t318 - mrSges(6,2) * t317) * t133 + t315 * t185 + t66 * t23 + t67 * t24 + (Ifges(7,5) * t88 + Ifges(7,6) * t89) * t364 + (Ifges(7,5) * t117 + Ifges(7,6) * t116) * t365 + (Ifges(7,1) * t88 + Ifges(7,4) * t89) * t377 + (Ifges(7,1) * t117 + Ifges(7,4) * t116) * t378 + (Ifges(7,4) * t88 + Ifges(7,2) * t89) * t379 + (Ifges(7,4) * t117 + Ifges(7,2) * t116) * t380 + (Ifges(7,4) * t156 + Ifges(7,2) * t155) * t381 + (Ifges(7,1) * t156 + Ifges(7,4) * t155) * t382; -t289 * t123 + t154 * t124 - t218 * t189 + t219 * t190 - t402 * t70 + t87 * t71 + (t13 * t87 - t14 * t402 + t60) * m(7) + (t154 * t44 - t289 * t45 + t118) * m(6) + (t109 * t219 - t110 * t218 + t143) * m(5) + t384; (-t154 * t46 + t264 * t23 + t260 * t24 + (-t260 * t71 + t264 * t70) * qJD(6) + (-t154 * t79 + t2 * t260 + t264 * t3 + (-t13 * t260 + t14 * t264) * qJD(6)) * m(7)) * pkin(5) + t300 + (Ifges(6,5) * t289 - Ifges(6,6) * t154) * t363 - t133 * (mrSges(6,1) * t154 + mrSges(6,2) * t289) + (Ifges(6,1) * t289 - t347) * t369 + (t154 * t45 + t289 * t44) * mrSges(6,3) + (-Ifges(6,2) * t154 + t150 + t75) * t371 - m(7) * (t13 * t15 + t14 * t16) + t87 * t400 - t385 - t44 * t123 + t45 * t124 - t16 * t70 - t15 * t71 + t74 * t368 + t410; -t13 * t70 + t14 * t71 + t38 * t377 + t410;];
tauc  = t1(:);

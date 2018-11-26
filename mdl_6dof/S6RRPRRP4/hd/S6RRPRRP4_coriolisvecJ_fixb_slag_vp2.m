% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRP4
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:12:28
% EndTime: 2018-11-23 17:12:46
% DurationCPUTime: 17.47s
% Computational Cost: add. (12703->661), mult. (32230->867), div. (0->0), fcn. (23539->8), ass. (0->301)
t237 = sin(qJ(4));
t290 = qJD(4) * t237;
t238 = sin(qJ(2));
t241 = cos(qJ(2));
t310 = sin(pkin(10));
t266 = qJD(1) * t310;
t311 = cos(pkin(10));
t267 = qJD(1) * t311;
t199 = -t238 * t266 + t241 * t267;
t305 = t199 * t237;
t424 = t290 - t305;
t410 = Ifges(6,1) + Ifges(7,1);
t409 = -Ifges(6,4) + Ifges(7,5);
t408 = Ifges(7,4) + Ifges(6,5);
t407 = -Ifges(6,6) + Ifges(7,6);
t406 = -Ifges(6,3) - Ifges(7,2);
t239 = cos(qJ(5));
t236 = sin(qJ(5));
t240 = cos(qJ(4));
t298 = t236 * t240;
t215 = t237 * t239 + t298;
t138 = t215 * t199;
t385 = qJD(4) + qJD(5);
t166 = t385 * t215;
t391 = t138 - t166;
t299 = t236 * t237;
t214 = -t239 * t240 + t299;
t139 = t214 * t199;
t165 = t385 * t214;
t390 = t139 - t165;
t194 = qJD(4) - t199;
t187 = qJD(5) + t194;
t348 = t187 / 0.2e1;
t200 = -t238 * t267 - t241 * t266;
t172 = qJD(2) * t237 - t200 * t240;
t247 = qJD(2) * t240 + t200 * t237;
t243 = t239 * t172 + t236 * t247;
t357 = t243 / 0.2e1;
t112 = t236 * t172 - t239 * t247;
t360 = t112 / 0.2e1;
t361 = -t112 / 0.2e1;
t330 = -qJ(3) - pkin(7);
t223 = t330 * t241;
t217 = qJD(1) * t223;
t203 = t310 * t217;
t222 = t330 * t238;
t216 = qJD(1) * t222;
t208 = qJD(2) * pkin(2) + t216;
t160 = t208 * t311 + t203;
t153 = -qJD(2) * pkin(3) - t160;
t108 = -pkin(4) * t247 + t153;
t47 = t112 * pkin(5) - qJ(6) * t243 + t108;
t413 = mrSges(7,3) * t47;
t423 = Ifges(6,4) * t361 + Ifges(7,5) * t360 + t348 * t408 + t357 * t410 - t413;
t110 = Ifges(6,4) * t112;
t316 = Ifges(7,5) * t112;
t402 = t408 * t187 + t243 * t410 - t110 + t316;
t422 = t402 / 0.2e1;
t278 = t310 * pkin(2);
t229 = t278 + pkin(8);
t331 = pkin(9) + t229;
t270 = qJD(4) * t331;
t196 = t237 * t270;
t209 = t331 * t237;
t210 = t331 * t240;
t248 = -t239 * t209 - t210 * t236;
t263 = t240 * t270;
t304 = t199 * t240;
t293 = qJD(1) * t238;
t283 = pkin(2) * t293;
t143 = -pkin(3) * t200 - pkin(8) * t199 + t283;
t163 = t216 * t311 + t203;
t92 = t240 * t143 - t163 * t237;
t74 = -pkin(4) * t200 - pkin(9) * t304 + t92;
t93 = t237 * t143 + t240 * t163;
t86 = -pkin(9) * t305 + t93;
t399 = qJD(5) * t248 + (-t196 - t86) * t239 + (-t263 - t74) * t236;
t212 = t238 * t311 + t241 * t310;
t289 = qJD(4) * t240;
t211 = t238 * t310 - t241 * t311;
t202 = t211 * qJD(2);
t303 = t202 * t237;
t421 = t212 * t289 - t303;
t269 = t311 * t217;
t162 = t216 * t310 - t269;
t388 = pkin(4) * t424 - t162;
t234 = -pkin(2) * t241 - pkin(1);
t294 = qJD(1) * t234;
t218 = qJD(3) + t294;
t132 = -t199 * pkin(3) + t200 * pkin(8) + t218;
t161 = t310 * t208 - t269;
t154 = qJD(2) * pkin(8) + t161;
t88 = t237 * t132 + t240 * t154;
t76 = pkin(9) * t247 + t88;
t313 = t239 * t76;
t87 = t240 * t132 - t154 * t237;
t75 = -pkin(9) * t172 + t87;
t66 = pkin(4) * t194 + t75;
t24 = t236 * t66 + t313;
t20 = qJ(6) * t187 + t24;
t420 = mrSges(6,1) * t108 - mrSges(7,2) * t20 - mrSges(6,3) * t24;
t77 = pkin(5) * t243 + qJ(6) * t112;
t109 = Ifges(7,5) * t243;
t60 = Ifges(7,6) * t187 + Ifges(7,3) * t112 + t109;
t317 = Ifges(6,4) * t243;
t63 = -Ifges(6,2) * t112 + Ifges(6,6) * t187 + t317;
t419 = -t47 * mrSges(7,1) - t60 / 0.2e1 + t63 / 0.2e1;
t315 = t236 * t76;
t23 = t239 * t66 - t315;
t393 = qJD(6) - t23;
t19 = -pkin(5) * t187 + t393;
t418 = mrSges(6,2) * t108 + mrSges(7,2) * t19 - mrSges(6,3) * t23 + t422;
t394 = Ifges(4,6) * qJD(2);
t396 = t247 * Ifges(5,6);
t417 = -t24 * mrSges(6,2) + t20 * mrSges(7,3) + t218 * mrSges(4,1) + t23 * mrSges(6,1) + t87 * mrSges(5,1) - t394 / 0.2e1 - t19 * mrSges(7,1) + t396 / 0.2e1 - t88 * mrSges(5,2);
t189 = qJD(1) * t202;
t245 = t247 * qJD(4);
t126 = -t189 * t240 + t245;
t127 = -qJD(4) * t172 + t189 * t237;
t53 = -qJD(5) * t112 + t239 * t126 + t236 * t127;
t369 = t53 / 0.2e1;
t54 = qJD(5) * t243 + t236 * t126 - t239 * t127;
t367 = t54 / 0.2e1;
t201 = t212 * qJD(2);
t188 = qJD(1) * t201;
t347 = t188 / 0.2e1;
t416 = -t247 / 0.2e1;
t405 = t408 * t188 + t409 * t54 + t410 * t53;
t43 = mrSges(6,1) * t188 - mrSges(6,3) * t53;
t44 = -t188 * mrSges(7,1) + t53 * mrSges(7,2);
t404 = t44 - t43;
t45 = -mrSges(6,2) * t188 - mrSges(6,3) * t54;
t46 = -mrSges(7,2) * t54 + mrSges(7,3) * t188;
t403 = t45 + t46;
t401 = -pkin(5) * t391 - qJ(6) * t390 - qJD(6) * t215 + t388;
t323 = mrSges(6,3) * t243;
t97 = mrSges(6,1) * t187 - t323;
t98 = -mrSges(7,1) * t187 + mrSges(7,2) * t243;
t327 = t98 - t97;
t400 = qJ(6) * t200 + t399;
t397 = Ifges(5,3) * t194;
t395 = Ifges(4,5) * qJD(2);
t392 = Ifges(5,5) * t126 + Ifges(5,6) * t127;
t159 = t211 * pkin(3) - t212 * pkin(8) + t234;
t170 = t222 * t310 - t223 * t311;
t164 = t240 * t170;
t107 = t237 * t159 + t164;
t325 = mrSges(4,3) * t200;
t295 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t247 + t172 * mrSges(5,2) - t325;
t171 = Ifges(5,4) * t247;
t103 = t172 * Ifges(5,1) + t194 * Ifges(5,5) + t171;
t296 = t240 * t103;
t320 = Ifges(5,4) * t172;
t102 = Ifges(5,2) * t247 + Ifges(5,6) * t194 + t320;
t297 = t237 * t102;
t389 = t296 / 0.2e1 - t297 / 0.2e1;
t387 = -t188 * t406 + t407 * t54 + t408 * t53;
t271 = qJD(2) * t330;
t197 = qJD(3) * t241 + t238 * t271;
t176 = t197 * qJD(1);
t198 = -t238 * qJD(3) + t241 * t271;
t244 = qJD(1) * t198;
t125 = t176 * t311 + t244 * t310;
t286 = qJD(1) * qJD(2);
t274 = t238 * t286;
t264 = pkin(2) * t274;
t129 = pkin(3) * t188 + pkin(8) * t189 + t264;
t39 = t240 * t125 + t237 * t129 + t132 * t289 - t154 * t290;
t40 = -t88 * qJD(4) - t125 * t237 + t240 * t129;
t386 = -t237 * t40 + t240 * t39;
t292 = qJD(1) * t241;
t308 = Ifges(3,6) * qJD(2);
t322 = Ifges(3,4) * t238;
t384 = t308 / 0.2e1 + (t241 * Ifges(3,2) + t322) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t292);
t383 = Ifges(5,5) * t172 + t112 * t407 - t187 * t406 + t243 * t408 + t396 + t397;
t22 = pkin(4) * t188 - pkin(9) * t126 + t40;
t26 = pkin(9) * t127 + t39;
t6 = -qJD(5) * t24 + t22 * t239 - t236 * t26;
t106 = t240 * t159 - t170 * t237;
t300 = t212 * t240;
t84 = pkin(4) * t211 - pkin(9) * t300 + t106;
t301 = t212 * t237;
t90 = -pkin(9) * t301 + t107;
t329 = t236 * t84 + t239 * t90;
t142 = t197 * t311 + t198 * t310;
t291 = qJD(2) * t238;
t144 = pkin(2) * t291 + pkin(3) * t201 + pkin(8) * t202;
t265 = -t142 * t237 + t240 * t144;
t302 = t202 * t240;
t34 = pkin(9) * t302 + pkin(4) * t201 + (-t164 + (pkin(9) * t212 - t159) * t237) * qJD(4) + t265;
t55 = t240 * t142 + t237 * t144 + t159 * t289 - t170 * t290;
t38 = -pkin(9) * t421 + t55;
t10 = -qJD(5) * t329 - t236 * t38 + t239 * t34;
t382 = t212 * t385;
t380 = -t395 / 0.2e1 - t218 * mrSges(4,2);
t349 = -t187 / 0.2e1;
t358 = -t243 / 0.2e1;
t378 = Ifges(6,4) * t360 + Ifges(7,5) * t361 + t349 * t408 + t358 * t410 + t413;
t124 = t176 * t310 - t311 * t244;
t85 = -t127 * pkin(4) + t124;
t11 = t54 * pkin(5) - t53 * qJ(6) - qJD(6) * t243 + t85;
t368 = -t54 / 0.2e1;
t377 = mrSges(6,1) * t85 + mrSges(7,1) * t11 + 0.2e1 * Ifges(7,3) * t367 - t53 * Ifges(6,4) / 0.2e1 - t188 * Ifges(6,6) / 0.2e1 + Ifges(7,6) * t347 + (t409 + Ifges(7,5)) * t369 + (-t368 + t367) * Ifges(6,2);
t376 = -Ifges(6,2) * t360 + Ifges(7,3) * t361 + t349 * t407 + t358 * t409 + t419;
t375 = -Ifges(6,2) * t361 + Ifges(7,3) * t360 + t348 * t407 + t357 * t409 - t419;
t373 = -0.2e1 * pkin(1);
t356 = t126 / 0.2e1;
t355 = t127 / 0.2e1;
t351 = -t172 / 0.2e1;
t350 = t172 / 0.2e1;
t346 = -t194 / 0.2e1;
t345 = -t199 / 0.2e1;
t344 = t199 / 0.2e1;
t343 = -t200 / 0.2e1;
t342 = t200 / 0.2e1;
t338 = -t237 / 0.2e1;
t337 = t240 / 0.2e1;
t336 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t293);
t334 = t39 * mrSges(5,2);
t333 = t40 * mrSges(5,1);
t95 = -mrSges(7,2) * t112 + mrSges(7,3) * t187;
t324 = mrSges(6,3) * t112;
t96 = -mrSges(6,2) * t187 - t324;
t328 = t95 + t96;
t326 = mrSges(4,3) * t199;
t321 = Ifges(4,4) * t200;
t319 = Ifges(5,4) * t237;
t318 = Ifges(5,4) * t240;
t309 = Ifges(3,5) * qJD(2);
t169 = -t311 * t222 - t223 * t310;
t307 = t124 * t169;
t288 = qJD(5) * t236;
t287 = qJD(5) * t239;
t280 = Ifges(5,3) * t188 + t392;
t279 = t311 * pkin(2);
t273 = t241 * t286;
t268 = t188 * mrSges(4,1) - t189 * mrSges(4,2);
t231 = -t279 - pkin(3);
t262 = mrSges(5,1) * t240 - mrSges(5,2) * t237;
t261 = mrSges(5,1) * t237 + mrSges(5,2) * t240;
t260 = Ifges(5,1) * t240 - t319;
t259 = Ifges(5,1) * t237 + t318;
t258 = -Ifges(5,2) * t237 + t318;
t257 = Ifges(5,2) * t240 + t319;
t256 = Ifges(5,5) * t240 - Ifges(5,6) * t237;
t255 = Ifges(5,5) * t237 + Ifges(5,6) * t240;
t31 = -t236 * t86 + t239 * t74;
t41 = -t236 * t90 + t239 * t84;
t251 = -t237 * t39 - t240 * t40;
t250 = t237 * t87 - t240 * t88;
t141 = t197 * t310 - t311 * t198;
t130 = -t194 * mrSges(5,2) + mrSges(5,3) * t247;
t131 = mrSges(5,1) * t194 - mrSges(5,3) * t172;
t249 = t130 * t240 - t131 * t237;
t157 = -t209 * t236 + t210 * t239;
t140 = pkin(4) * t301 + t169;
t5 = t236 * t22 + t239 * t26 + t66 * t287 - t288 * t76;
t9 = t236 * t34 + t239 * t38 + t84 * t287 - t288 * t90;
t246 = t153 * t261;
t219 = -t240 * pkin(4) + t231;
t94 = pkin(4) * t421 + t141;
t2 = qJ(6) * t188 + qJD(6) * t187 + t5;
t3 = -pkin(5) * t188 - t6;
t242 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t387;
t235 = Ifges(3,4) * t292;
t233 = -pkin(4) * t239 - pkin(5);
t230 = pkin(4) * t236 + qJ(6);
t224 = pkin(4) * t287 + qJD(6);
t207 = Ifges(3,1) * t293 + t235 + t309;
t192 = Ifges(4,4) * t199;
t174 = -qJD(2) * mrSges(4,2) + t326;
t155 = -mrSges(4,1) * t199 - mrSges(4,2) * t200;
t150 = t214 * pkin(5) - t215 * qJ(6) + t219;
t149 = t214 * t212;
t148 = t215 * t212;
t147 = -t200 * Ifges(4,1) + t192 + t395;
t146 = t199 * Ifges(4,2) - t321 + t394;
t105 = -mrSges(5,2) * t188 + mrSges(5,3) * t127;
t104 = mrSges(5,1) * t188 - mrSges(5,3) * t126;
t100 = qJD(5) * t157 - t196 * t236 + t239 * t263;
t83 = -mrSges(5,1) * t127 + mrSges(5,2) * t126;
t79 = mrSges(6,1) * t112 + mrSges(6,2) * t243;
t78 = mrSges(7,1) * t112 - mrSges(7,3) * t243;
t71 = t126 * Ifges(5,1) + t127 * Ifges(5,4) + t188 * Ifges(5,5);
t70 = t126 * Ifges(5,4) + t127 * Ifges(5,2) + t188 * Ifges(5,6);
t69 = t148 * pkin(5) + t149 * qJ(6) + t140;
t68 = -t202 * t298 + (t300 * t385 - t303) * t239 - t299 * t382;
t67 = t214 * t202 - t215 * t382;
t58 = pkin(4) * t172 + t77;
t56 = -qJD(4) * t107 + t265;
t36 = -pkin(5) * t211 - t41;
t35 = qJ(6) * t211 + t329;
t30 = t239 * t75 - t315;
t29 = t236 * t75 + t313;
t28 = pkin(5) * t200 - t31;
t18 = mrSges(6,1) * t54 + mrSges(6,2) * t53;
t17 = mrSges(7,1) * t54 - mrSges(7,3) * t53;
t12 = t68 * pkin(5) - t67 * qJ(6) + t149 * qJD(6) + t94;
t8 = -pkin(5) * t201 - t10;
t7 = qJ(6) * t201 + qJD(6) * t211 + t9;
t1 = [(t280 + t387 + t392) * t211 / 0.2e1 - (t147 / 0.2e1 - t160 * mrSges(4,3) + Ifges(4,1) * t343 + Ifges(4,4) * t344 - t380 + t389) * t202 + (-t308 / 0.2e1 + (mrSges(3,1) * t373 - 0.3e1 / 0.2e1 * t322 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t241) * qJD(1) + (qJD(1) * (mrSges(4,1) * t211 + mrSges(4,2) * t212) + m(4) * (t218 + t294) + t155) * pkin(2) - t384) * t291 + (t87 * t302 + t88 * t303 + (qJD(4) * t250 + t251) * t212) * mrSges(5,3) + t247 * (-Ifges(5,4) * t302 + Ifges(5,2) * t303) / 0.2e1 + (-Ifges(5,1) * t302 + Ifges(5,4) * t303) * t350 + t194 * (-Ifges(5,5) * t302 + Ifges(5,6) * t303) / 0.2e1 + t153 * (-mrSges(5,1) * t303 - mrSges(5,2) * t302) + (-Ifges(7,5) * t149 + Ifges(7,6) * t211) * t367 + (-Ifges(6,4) * t149 + Ifges(6,6) * t211) * t368 + t6 * (mrSges(6,1) * t211 + mrSges(6,3) * t149) + t3 * (-mrSges(7,1) * t211 - mrSges(7,2) * t149) + m(6) * (t10 * t23 + t108 * t94 + t140 * t85 + t24 * t9 + t329 * t5 + t41 * t6) + t329 * t45 + t234 * t268 + (-t125 * t211 - t170 * t188) * mrSges(4,3) + (t71 * t337 + t70 * t338 - Ifges(4,1) * t189 - Ifges(4,4) * t188 + t258 * t355 + t260 * t356 + (mrSges(4,3) + t261) * t124 + (t103 * t338 - t240 * t102 / 0.2e1 + t257 * t416 + t255 * t346 + t259 * t351 + t153 * t262) * qJD(4)) * t212 - (mrSges(4,3) * t169 - Ifges(4,4) * t211) * t189 + (t418 + t423) * t67 - t405 * t149 / 0.2e1 + (-mrSges(7,2) * t2 - mrSges(6,3) * t5 + t347 * t407 + t377) * t148 + (-t149 * t408 - t211 * t406 + t212 * t256) * t347 + (-t149 * t410 + t211 * t408) * t369 + (Ifges(4,2) + Ifges(5,3) / 0.2e1) * t188 * t211 + m(7) * (t11 * t69 + t12 * t47 + t19 * t8 + t2 * t35 + t20 * t7 + t3 * t36) + (t207 / 0.2e1 - t336 + t309 / 0.2e1 + (mrSges(3,2) * t373 + 0.3e1 / 0.2e1 * Ifges(3,4) * t241) * qJD(1)) * qJD(2) * t241 + t142 * t174 + t169 * t83 + t140 * t18 + t55 * t130 + t56 * t131 + t106 * t104 + t107 * t105 + t8 * t98 + t94 * t79 + t7 * t95 + t9 * t96 + t10 * t97 + t12 * t78 + t69 * t17 + t41 * t43 + t36 * t44 + t35 * t46 + (t375 + t420) * t68 + (t383 / 0.2e1 - t406 * t348 + t408 * t357 - t161 * mrSges(4,3) - Ifges(4,4) * t343 - Ifges(4,2) * t344 + Ifges(5,5) * t350 + Ifges(7,6) * t360 + Ifges(6,6) * t361 + t397 / 0.2e1 - t146 / 0.2e1 + t417) * t201 + t295 * t141 + t211 * t333 + m(4) * (t125 * t170 - t141 * t160 + t142 * t161 + t307) + m(5) * (t106 * t40 + t107 * t39 + t141 * t153 + t55 * t88 + t56 * t87 + t307) + (t11 * t149 + t2 * t211) * mrSges(7,3) - t211 * t334 + (-t149 * t85 - t211 * t5) * mrSges(6,2); t388 * t79 + (t246 + t389) * qJD(4) + (-m(5) * t153 - t295) * t162 + (-mrSges(6,1) * t391 + mrSges(6,2) * t390) * t108 + (-t131 * t289 - t130 * t290 + t240 * t105 - t237 * t104 + m(5) * ((-t237 * t88 - t240 * t87) * qJD(4) + t386)) * t229 + (t321 + t383) * t342 + t384 * t293 + t375 * t166 + t376 * t138 + t377 * t214 - (-Ifges(3,2) * t293 + t207 + t235) * t292 / 0.2e1 + (t172 * t260 + t194 * t256) * qJD(4) / 0.2e1 + (pkin(1) * (mrSges(3,1) * t238 + mrSges(3,2) * t241) - t238 * (Ifges(3,1) * t241 - t322) / 0.2e1) * qJD(1) ^ 2 + (t192 + t296 + t147) * t345 + (t248 * t6 + t157 * t5 + t219 * t85 + t399 * t24 + (-t100 - t31) * t23 + t388 * t108) * m(6) + (t11 * t150 - t248 * t3 + t157 * t2 + t401 * t47 + t400 * t20 + (t100 - t28) * t19) * m(7) - t404 * t248 + (-t188 * t278 + t189 * t279) * mrSges(4,3) + Ifges(3,5) * t273 + (Ifges(4,1) * t342 + t256 * t346 + t258 * t416 + t260 * t351 - t246 + t380) * t199 - (t422 + t423) * t165 + (-mrSges(3,1) * t273 + mrSges(3,2) * t274) * pkin(7) + t399 * t96 + t400 * t95 + t327 * t100 + t401 * t78 + t403 * t157 + t258 * t245 / 0.2e1 - m(5) * (t87 * t92 + t88 * t93) + (m(5) * t231 - mrSges(4,1) - t262) * t124 - Ifges(3,6) * t274 + (-t424 * t88 + (-t289 + t304) * t87 + t386) * mrSges(5,3) + ((-t124 * t311 + t125 * t310) * pkin(2) + t160 * t162 - t161 * t163 - t218 * t283) * m(4) - Ifges(4,5) * t189 - Ifges(4,6) * t188 - t163 * t174 + t150 * t17 - t93 * t130 - t92 * t131 - t125 * mrSges(4,2) - t31 * t97 - t28 * t98 + (t214 * t407 + t255) * t347 + (t3 * mrSges(7,2) - t6 * mrSges(6,3) + t405 / 0.2e1 + t408 * t347 + t85 * mrSges(6,2) - t11 * mrSges(7,3) + Ifges(6,4) * t368 + Ifges(7,5) * t367 + t410 * t369) * t215 + (t19 * t390 - t2 * t214 + t20 * t391) * mrSges(7,2) + (-t214 * t5 - t23 * t390 + t24 * t391) * mrSges(6,3) + (-t378 + t422) * t139 + (-Ifges(5,5) * t351 + Ifges(4,2) * t345 - Ifges(6,6) * t360 - Ifges(7,6) * t361 - Ifges(5,3) * t346 + t349 * t406 - t358 * t408 + t417) * t200 - t155 * t283 + t219 * t18 + t231 * t83 - (Ifges(3,5) * t241 - Ifges(3,6) * t238) * t286 / 0.2e1 + t237 * t71 / 0.2e1 + t160 * t326 + t292 * t336 + t70 * t337 + t146 * t343 + t297 * t344 + t257 * t355 + t259 * t356 - t161 * t325; t240 * t104 + t237 * t105 + t403 * t215 + t404 * t214 + t249 * qJD(4) + (-t174 - t249) * t199 + (t78 + t79 + t295) * t200 + t268 + t390 * t328 - t391 * t327 + (-t19 * t391 + t2 * t215 + t20 * t390 + t200 * t47 + t214 * t3) * m(7) + (t108 * t200 - t214 * t6 + t215 * t5 + t23 * t391 + t24 * t390) * m(6) + (t153 * t200 - t194 * t250 - t251) * m(5) + (-t160 * t200 - t161 * t199 + t264) * m(4); (t376 - t420) * t243 - t153 * (t172 * mrSges(5,1) + mrSges(5,2) * t247) + (t88 * t172 + t247 * t87) * mrSges(5,3) + t242 + t333 + t280 - t87 * t130 + t88 * t131 - t58 * t78 + t224 * t95 + t230 * t46 + t233 * t44 - t334 + (Ifges(5,5) * t247 - Ifges(5,6) * t172) * t346 + t102 * t350 + (Ifges(5,1) * t247 - t320) * t351 - m(6) * (-t23 * t29 + t24 * t30) + (-t19 * t29 + t2 * t230 + t233 * t3 - t47 * t58 + (t224 - t30) * t20) * m(7) - t327 * t29 - t328 * t30 + (-t172 * t79 + t236 * t45 + t239 * t43 + (t236 * t327 + t239 * t96) * qJD(5) + m(7) * t19 * t288 + (0.2e1 * t108 * t351 - t23 * t288 + t236 * t5 + t239 * t6 + t24 * t287) * m(6)) * pkin(4) + (-Ifges(5,2) * t172 + t103 + t171) * t416 + (-t378 + t418) * t112; t242 + (t323 - t327) * t24 + (-t324 - t328) * t23 + (t112 * t19 + t20 * t243) * mrSges(7,2) - t108 * (mrSges(6,1) * t243 - mrSges(6,2) * t112) + t63 * t357 + (Ifges(7,3) * t243 - t316) * t361 - t47 * (mrSges(7,1) * t243 + mrSges(7,3) * t112) + qJD(6) * t95 - t77 * t78 - pkin(5) * t44 + qJ(6) * t46 + (-t112 * t408 + t243 * t407) * t349 + (-pkin(5) * t3 + qJ(6) * t2 - t19 * t24 + t20 * t393 - t47 * t77) * m(7) + (-Ifges(6,2) * t243 - t110 + t402) * t360 + (-t112 * t410 + t109 - t317 + t60) * t358; t243 * t78 - t187 * t95 + 0.2e1 * (t3 / 0.2e1 + t47 * t357 + t20 * t349) * m(7) + t44;];
tauc  = t1(:);

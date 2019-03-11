% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:27
% EndTime: 2019-03-09 12:45:51
% DurationCPUTime: 12.78s
% Computational Cost: add. (8130->594), mult. (18639->787), div. (0->0), fcn. (11496->6), ass. (0->271)
t240 = sin(qJ(4));
t243 = cos(qJ(4));
t244 = cos(qJ(2));
t305 = qJD(1) * t244;
t186 = -qJD(2) * t240 - t243 * t305;
t303 = qJD(2) * t243;
t187 = -t240 * t305 + t303;
t239 = sin(qJ(5));
t242 = cos(qJ(5));
t123 = t186 * t239 + t187 * t242;
t229 = pkin(7) * t305;
t195 = pkin(3) * t305 + t229;
t238 = qJD(2) * qJ(3);
t171 = t238 + t195;
t130 = -pkin(4) * t186 + t171;
t241 = sin(qJ(2));
t233 = t241 * qJD(1);
t220 = t233 + qJD(4);
t245 = -pkin(2) - pkin(8);
t278 = -qJ(3) * t241 - pkin(1);
t178 = t244 * t245 + t278;
t147 = t178 * qJD(1);
t227 = pkin(7) * t233;
t194 = -pkin(3) * t233 - t227;
t365 = qJD(3) - t194;
t152 = qJD(2) * t245 + t365;
t92 = -t147 * t240 + t152 * t243;
t80 = -pkin(9) * t187 + t92;
t71 = pkin(4) * t220 + t80;
t93 = t147 * t243 + t152 * t240;
t81 = pkin(9) * t186 + t93;
t77 = t239 * t81;
t22 = t242 * t71 - t77;
t405 = qJ(6) * t123;
t18 = t22 - t405;
t293 = qJD(4) + qJD(5);
t211 = t233 + t293;
t17 = pkin(5) * t211 + t18;
t79 = t242 * t81;
t23 = t239 * t71 + t79;
t275 = t186 * t242 - t187 * t239;
t369 = qJ(6) * t275;
t19 = t23 + t369;
t295 = qJD(1) * qJD(2);
t279 = t244 * t295;
t216 = Ifges(7,3) * t279;
t217 = Ifges(6,3) * t279;
t339 = -t211 / 0.2e1;
t346 = t123 / 0.2e1;
t347 = -t123 / 0.2e1;
t294 = qJD(2) * qJD(4);
t298 = qJD(4) * t244;
t304 = qJD(2) * t241;
t138 = -t240 * t294 + (t240 * t304 - t243 * t298) * qJD(1);
t283 = t240 * t298;
t250 = t241 * t303 + t283;
t139 = qJD(1) * t250 - t243 * t294;
t56 = qJD(5) * t275 + t138 * t242 + t139 * t239;
t280 = t241 * t295;
t219 = pkin(2) * t280;
t262 = pkin(8) * t241 - qJ(3) * t244;
t301 = qJD(3) * t241;
t249 = qJD(2) * t262 - t301;
t134 = qJD(1) * t249 + t219;
t302 = qJD(2) * t244;
t355 = pkin(3) + pkin(7);
t197 = t355 * t302;
t177 = qJD(1) * t197;
t47 = -qJD(4) * t93 - t134 * t240 + t177 * t243;
t27 = pkin(4) * t279 - pkin(9) * t138 + t47;
t299 = qJD(4) * t243;
t300 = qJD(4) * t240;
t46 = t134 * t243 - t147 * t300 + t152 * t299 + t177 * t240;
t33 = pkin(9) * t139 + t46;
t6 = -qJD(5) * t23 - t239 * t33 + t242 * t27;
t2 = pkin(5) * t279 - qJ(6) * t56 - qJD(6) * t123 + t6;
t296 = qJD(5) * t242;
t297 = qJD(5) * t239;
t5 = t239 * t27 + t242 * t33 + t296 * t71 - t297 * t81;
t57 = -qJD(5) * t123 - t138 * t239 + t139 * t242;
t3 = qJ(6) * t57 + qJD(6) * t275 + t5;
t364 = mrSges(6,1) * t6 + mrSges(7,1) * t2 - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t392 = Ifges(6,5) + Ifges(7,5);
t394 = Ifges(6,1) + Ifges(7,1);
t393 = Ifges(6,4) + Ifges(7,4);
t413 = t393 * t275;
t384 = t123 * t394 + t211 * t392 + t413;
t390 = Ifges(6,6) + Ifges(7,6);
t391 = Ifges(6,2) + Ifges(7,2);
t410 = t123 * t393;
t385 = t211 * t390 + t275 * t391 + t410;
t398 = -t275 / 0.2e1;
t51 = Ifges(7,6) * t57;
t52 = Ifges(6,6) * t57;
t53 = Ifges(7,5) * t56;
t54 = Ifges(6,5) * t56;
t74 = -pkin(5) * t275 + qJD(6) + t130;
t416 = t216 + t217 + t51 + t52 + t53 + t54 + t364 + (-t123 * t390 + t275 * t392) * t339 + (t123 * t19 + t17 * t275) * mrSges(7,3) + (t123 * t23 + t22 * t275) * mrSges(6,3) - t130 * (mrSges(6,1) * t123 + mrSges(6,2) * t275) - t74 * (mrSges(7,1) * t123 + mrSges(7,2) * t275) + t385 * t346 + (-t123 * t391 + t384 + t413) * t398 + (t275 * t394 - t410) * t347;
t228 = pkin(2) * t233;
t156 = qJD(1) * t262 + t228;
t113 = -t156 * t240 + t195 * t243;
t255 = -pkin(9) * t240 * t241 + pkin(4) * t244;
t332 = pkin(9) - t245;
t415 = -qJD(1) * t255 + t300 * t332 - t113;
t114 = t156 * t243 + t195 * t240;
t201 = t332 * t243;
t284 = t243 * t233;
t414 = pkin(9) * t284 + qJD(4) * t201 + t114;
t310 = t242 * t243;
t127 = -t239 * t300 - t240 * t297 + t293 * t310;
t311 = t239 * t240;
t148 = -t233 * t311 + t242 * t284;
t308 = t127 + t148;
t257 = t239 * t243 + t240 * t242;
t128 = t293 * t257;
t251 = t257 * t241;
t149 = qJD(1) * t251;
t412 = t128 + t149;
t318 = t187 * Ifges(5,4);
t110 = Ifges(5,2) * t186 + Ifges(5,6) * t220 + t318;
t179 = Ifges(5,4) * t186;
t111 = Ifges(5,1) * t187 + Ifges(5,5) * t220 + t179;
t260 = t240 * t92 - t243 * t93;
t328 = Ifges(5,4) * t240;
t264 = Ifges(5,2) * t243 + t328;
t327 = Ifges(5,4) * t243;
t266 = Ifges(5,1) * t240 + t327;
t269 = mrSges(5,1) * t243 - mrSges(5,2) * t240;
t320 = Ifges(5,6) * t243;
t324 = Ifges(5,5) * t240;
t335 = -t243 / 0.2e1;
t336 = -t240 / 0.2e1;
t337 = -t220 / 0.2e1;
t342 = -t187 / 0.2e1;
t343 = -t186 / 0.2e1;
t411 = t260 * mrSges(5,3) + t110 * t335 + t111 * t336 + t171 * t269 + (t320 + t324) * t337 + t264 * t343 + t266 * t342;
t409 = qJD(3) + t227;
t407 = -t393 * t57 / 0.2e1 - t394 * t56 / 0.2e1 - t392 * t279 / 0.2e1;
t406 = -mrSges(3,1) + mrSges(4,2);
t200 = t332 * t240;
t383 = t200 * t297 - t201 * t296 + t239 * t415 - t242 * t414;
t132 = -t200 * t242 - t201 * t239;
t382 = -qJD(5) * t132 + t239 * t414 + t242 * t415;
t285 = -pkin(4) * t243 - pkin(3);
t367 = pkin(4) * t299 - t233 * t285 + t409;
t202 = -pkin(2) * t244 + t278;
t172 = t202 * qJD(1);
t205 = -t229 - t238;
t395 = qJD(2) / 0.2e1;
t396 = -qJD(2) / 0.2e1;
t397 = -qJD(1) / 0.2e1;
t404 = Ifges(3,6) * t395 + (Ifges(3,4) * t241 + Ifges(3,2) * t244) * qJD(1) / 0.2e1 + Ifges(4,5) * t396 + (-Ifges(4,6) * t241 - Ifges(4,3) * t244) * t397 - t205 * mrSges(4,1) + t172 * mrSges(4,2) + t411;
t256 = -t310 + t311;
t402 = t22 * t412 - t23 * t308 + t256 * t6 - t257 * t5;
t401 = t17 * t412 - t19 * t308 + t2 * t256 - t257 * t3;
t338 = t211 / 0.2e1;
t349 = t275 / 0.2e1;
t389 = t279 * t390 + t391 * t57 + t393 * t56;
t387 = -qJ(6) * t308 - qJD(6) * t257 + t383;
t386 = -pkin(5) * t305 + qJ(6) * t412 + qJD(6) * t256 + t382;
t368 = pkin(5) * t308 + t367;
t209 = t355 * t241;
t191 = t243 * t209;
t277 = pkin(9) * t244 - t178;
t106 = pkin(4) * t241 + t240 * t277 + t191;
t190 = t240 * t209;
t126 = t178 * t243 + t190;
t309 = t243 * t244;
t112 = -pkin(9) * t309 + t126;
t60 = t106 * t239 + t112 * t242;
t366 = t240 * t46 + t243 * t47;
t363 = t244 * t293;
t362 = t47 * mrSges(5,1) - t46 * mrSges(5,2) + Ifges(5,5) * t138 + Ifges(5,6) * t139;
t199 = -qJD(2) * pkin(2) + t409;
t226 = Ifges(3,4) * t305;
t287 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t288 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t290 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t322 = Ifges(4,6) * t244;
t361 = t290 * t123 + t287 * t211 + t288 * t275 + t17 * mrSges(7,1) + t199 * mrSges(4,1) + t22 * mrSges(6,1) + t92 * mrSges(5,1) + t220 * Ifges(5,3) + t187 * Ifges(5,5) + t186 * Ifges(5,6) + Ifges(3,1) * t233 / 0.2e1 + Ifges(3,5) * t395 + t226 / 0.2e1 + Ifges(4,4) * t396 + (-t241 * Ifges(4,2) - t322) * t397 - t172 * mrSges(4,3) - t19 * mrSges(7,2) - t23 * mrSges(6,2) - t93 * mrSges(5,2) + t390 * t349 + t392 * t346 + (Ifges(7,3) + Ifges(6,3)) * t338;
t359 = t56 / 0.2e1;
t358 = t57 / 0.2e1;
t354 = pkin(1) * mrSges(3,1);
t353 = pkin(1) * mrSges(3,2);
t352 = t110 / 0.2e1;
t42 = -mrSges(7,2) * t279 + mrSges(7,3) * t57;
t43 = -mrSges(6,2) * t279 + mrSges(6,3) * t57;
t331 = t42 + t43;
t30 = t242 * t80 - t77;
t96 = -mrSges(7,2) * t211 + mrSges(7,3) * t275;
t97 = -mrSges(6,2) * t211 + mrSges(6,3) * t275;
t330 = t96 + t97;
t98 = mrSges(7,1) * t211 - mrSges(7,3) * t123;
t99 = mrSges(6,1) * t211 - mrSges(6,3) * t123;
t329 = t98 + t99;
t323 = Ifges(5,5) * t243;
t321 = Ifges(5,6) * t240;
t313 = qJD(2) * mrSges(3,2);
t222 = pkin(4) * t240 + qJ(3);
t129 = -mrSges(5,1) * t186 + mrSges(5,2) * t187;
t207 = -mrSges(4,1) * t305 - qJD(2) * mrSges(4,3);
t306 = -t207 + t129;
t210 = t355 * t244;
t292 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t291 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t289 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t286 = m(4) * pkin(7) + mrSges(4,1);
t165 = pkin(4) * t309 + t210;
t15 = -mrSges(7,1) * t57 + mrSges(7,2) * t56;
t29 = -t239 * t80 - t79;
t59 = t106 * t242 - t112 * t239;
t232 = pkin(2) * t304;
t143 = t232 + t249;
t276 = -t143 * t240 + t197 * t243;
t131 = t200 * t239 - t201 * t242;
t196 = t355 * t304;
t273 = m(4) * t199 + (mrSges(4,1) + mrSges(3,3)) * t233 + t406 * qJD(2);
t272 = m(4) * t205 - mrSges(3,3) * t305 + t207 + t313;
t268 = mrSges(5,1) * t240 + mrSges(5,2) * t243;
t267 = Ifges(5,1) * t243 - t328;
t265 = -Ifges(5,2) * t240 + t327;
t118 = mrSges(5,1) * t279 - mrSges(5,3) * t138;
t119 = -mrSges(5,2) * t279 + mrSges(5,3) * t139;
t259 = t118 * t243 + t119 * t240;
t141 = -mrSges(5,2) * t220 + mrSges(5,3) * t186;
t142 = mrSges(5,1) * t220 - mrSges(5,3) * t187;
t258 = t141 * t243 - t142 * t240;
t48 = t255 * qJD(2) + (t243 * t277 - t190) * qJD(4) + t276;
t67 = t143 * t243 - t178 * t300 + t197 * t240 + t209 * t299;
t55 = pkin(9) * t250 + t67;
t9 = t106 * t296 - t112 * t297 + t239 * t48 + t242 * t55;
t252 = -qJ(3) * t302 - t301;
t237 = qJD(2) * qJD(3);
t153 = -qJD(1) * t196 + t237;
t102 = -pkin(4) * t139 + t153;
t10 = -qJD(5) * t60 - t239 * t55 + t242 * t48;
t133 = -pkin(4) * t283 + (-pkin(7) + t285) * t304;
t225 = pkin(4) * t242 + pkin(5);
t218 = Ifges(5,3) * t279;
t198 = pkin(7) * t280 - t237;
t192 = (mrSges(4,2) * t244 - mrSges(4,3) * t241) * qJD(1);
t160 = t232 + t252;
t158 = t257 * t244;
t157 = t256 * t244;
t150 = pkin(5) * t257 + t222;
t145 = qJD(1) * t252 + t219;
t125 = -t178 * t240 + t191;
t115 = -pkin(5) * t157 + t165;
t104 = -qJ(6) * t257 + t132;
t103 = qJ(6) * t256 + t131;
t90 = pkin(4) * t187 + pkin(5) * t123;
t88 = -mrSges(5,1) * t139 + mrSges(5,2) * t138;
t85 = Ifges(5,1) * t138 + Ifges(5,4) * t139 + Ifges(5,5) * t279;
t84 = Ifges(5,4) * t138 + Ifges(5,2) * t139 + Ifges(5,6) * t279;
t83 = -t256 * t304 + t257 * t363;
t82 = qJD(2) * t251 + t256 * t363;
t73 = -mrSges(6,1) * t275 + mrSges(6,2) * t123;
t72 = -mrSges(7,1) * t275 + mrSges(7,2) * t123;
t68 = -qJD(4) * t126 + t276;
t58 = -pkin(5) * t83 + t133;
t41 = mrSges(6,1) * t279 - mrSges(6,3) * t56;
t40 = mrSges(7,1) * t279 - mrSges(7,3) * t56;
t37 = qJ(6) * t157 + t60;
t36 = pkin(5) * t241 + qJ(6) * t158 + t59;
t24 = -pkin(5) * t57 + t102;
t21 = t30 - t405;
t20 = t29 - t369;
t16 = -mrSges(6,1) * t57 + mrSges(6,2) * t56;
t8 = qJ(6) * t83 + qJD(6) * t157 + t9;
t7 = pkin(5) * t302 - qJ(6) * t82 + qJD(6) * t158 + t10;
t1 = [t210 * t88 + t160 * t192 - t196 * t129 + t165 * t16 + t67 * t141 + t68 * t142 + t133 * t73 + t130 * (-mrSges(6,1) * t83 + mrSges(6,2) * t82) + t125 * t118 + t126 * t119 + m(4) * (t145 * t202 + t160 * t172) + t158 * t407 + m(5) * (t125 * t47 + t126 * t46 + t153 * t210 - t171 * t196 + t67 * t93 + t68 * t92) + t102 * (-mrSges(6,1) * t157 - mrSges(6,2) * t158) + t24 * (-mrSges(7,1) * t157 - mrSges(7,2) * t158) + (t157 * t5 + t158 * t6 - t22 * t82 + t23 * t83) * mrSges(6,3) + (t157 * t3 + t158 * t2 - t17 * t82 + t19 * t83) * mrSges(7,3) + t115 * t15 + t7 * t98 + t10 * t99 + t8 * t96 + t9 * t97 + t74 * (-mrSges(7,1) * t83 + mrSges(7,2) * t82) + t58 * t72 + t59 * t41 + t60 * t43 + t37 * t42 + t36 * t40 + (t157 * t391 - t158 * t393) * t358 + (t391 * t83 + t393 * t82) * t349 + (t157 * t393 - t158 * t394) * t359 + (t393 * t83 + t394 * t82) * t346 + (t390 * t83 + t392 * t82) * t338 + t384 * t82 / 0.2e1 + t385 * t83 / 0.2e1 + t389 * t157 / 0.2e1 + (-t145 * mrSges(4,3) + t53 / 0.2e1 + t51 / 0.2e1 + t216 / 0.2e1 + t54 / 0.2e1 + t52 / 0.2e1 + t217 / 0.2e1 + t218 / 0.2e1 + t288 * t57 + t290 * t56 + t362 + t364) * t241 + (-t138 * t266 / 0.2e1 + t145 * mrSges(4,2) + t153 * t269 - t139 * t264 / 0.2e1 + t85 * t336 + t84 * t335 - t286 * t198 + (t240 * t47 - t243 * t46) * mrSges(5,3) + (t220 * (t321 - t323) / 0.2e1 + t265 * t343 + t267 * t342 - t171 * t268 + t240 * t352 + t111 * t335 + (t240 * t93 + t243 * t92) * mrSges(5,3)) * qJD(4) + (t273 * pkin(7) + ((-t324 / 0.2e1 - t320 / 0.2e1 - t289) * t244 - 0.2e1 * t353 - t202 * mrSges(4,3) - t290 * t158 + t288 * t157 + (0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + t286 * pkin(7) + t287) * t241) * qJD(1) + t292 * qJD(2) + t361) * qJD(2)) * t244 + (t291 * qJD(2) + t272 * pkin(7) + (-mrSges(4,2) * t202 + t241 * t289 - 0.2e1 * t354) * qJD(1) - t404) * t304 + m(7) * (t115 * t24 + t17 * t7 + t19 * t8 + t2 * t36 + t3 * t37 + t58 * t74) + m(6) * (t10 * t22 + t102 * t165 + t130 * t133 + t23 * t9 + t5 * t60 + t59 * t6); t222 * t16 - t194 * t129 - t198 * mrSges(4,3) + t150 * t15 - t114 * t141 - t113 * t142 + t132 * t43 + t131 * t41 + t411 * qJD(4) - t366 * mrSges(5,3) - t389 * t257 / 0.2e1 + t24 * (mrSges(7,1) * t257 - mrSges(7,2) * t256) + t102 * (mrSges(6,1) * t257 - mrSges(6,2) * t256) + (-t256 * t394 - t257 * t393) * t359 + (-t256 * t393 - t257 * t391) * t358 + (t148 * t391 + t149 * t393) * t398 + (-t127 * t393 - t128 * t394) * t346 + (-t127 * t390 - t128 * t392) * t338 + (-t127 * t391 - t128 * t393) * t349 + t384 * (-t149 / 0.2e1 - t128 / 0.2e1) + t401 * mrSges(7,3) + t402 * mrSges(6,3) + t256 * t407 + (-m(4) * t172 - t192) * (-qJ(3) * t305 + t228) + m(4) * (-qJ(3) * t198 - qJD(3) * t205) + t103 * t40 + t104 * t42 + qJ(3) * t88 + (mrSges(6,1) * t308 - mrSges(6,2) * t412) * t130 + (mrSges(7,1) * t308 - mrSges(7,2) * t412) * t74 + (t148 * t393 + t149 * t394) * t347 + (t148 * t390 + t149 * t392) * t339 + t382 * t99 + t383 * t97 + (t102 * t222 + t130 * t367 + t131 * t6 + t132 * t5 + t22 * t382 + t23 * t383) * m(6) + t385 * (-t148 / 0.2e1 - t127 / 0.2e1) + t386 * t98 + t387 * t96 + (t103 * t2 + t104 * t3 + t150 * t24 + t17 * t386 + t19 * t387 + t368 * t74) * m(7) + t368 * t72 + ((-m(5) * t260 + t258) * qJD(4) + m(5) * t366 + t259) * t245 + t367 * t73 + (t153 * qJ(3) - t113 * t92 - t114 * t93 + t171 * t365) * m(5) + t306 * qJD(3) + t243 * t85 / 0.2e1 + (((-mrSges(4,1) * qJ(3) + t291) * qJD(2) + (-t272 + t313) * pkin(7) + (t354 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t241) * qJD(1) + t404) * t241 + ((-t322 / 0.2e1 + t353) * qJD(1) + ((-m(4) * pkin(2) + t406) * qJD(2) - t273) * pkin(7) + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t233 + (t323 / 0.2e1 - t321 / 0.2e1 - pkin(2) * mrSges(4,1) - t290 * t256 - t288 * t257 + t292) * qJD(2) - t226 / 0.2e1 - t361) * t244) * qJD(1) + t84 * t336 + t139 * t265 / 0.2e1 + t138 * t267 / 0.2e1 + t153 * t268; -(t40 + t41) * t256 + t331 * t257 + t258 * qJD(4) + (t192 + t258) * t233 + (t286 * t305 - t306 - t72 - t73) * qJD(2) - m(4) * (-qJD(2) * t205 - t172 * t233) + t259 + t308 * t330 - t412 * t329 + (-qJD(2) * t74 - t401) * m(7) + (-qJD(2) * t130 - t402) * m(6) + (-qJD(2) * t171 - t220 * t260 + t366) * m(5); t225 * t40 - t171 * (mrSges(5,1) * t187 + mrSges(5,2) * t186) - t92 * t141 + t93 * t142 + (t186 * t92 + t187 * t93) * mrSges(5,3) - m(6) * (t22 * t29 + t23 * t30) + (-Ifges(5,2) * t187 + t111 + t179) * t343 + (-t187 * t73 + t242 * t41 + t331 * t239 + (-t239 * t329 + t242 * t330) * qJD(5) + (-t130 * t187 - t22 * t297 + t23 * t296 + t239 * t5 + t242 * t6) * m(6)) * pkin(4) + (t2 * t225 + (-t17 * t297 + t19 * t296 + t239 * t3) * pkin(4) - t17 * t20 - t19 * t21 - t74 * t90) * m(7) + t362 - t20 * t98 - t29 * t99 - t90 * t72 - t21 * t96 - t30 * t97 + t218 + (Ifges(5,5) * t186 - Ifges(5,6) * t187) * t337 + (Ifges(5,1) * t186 - t318) * t342 + t187 * t352 + t416; (-(-t17 + t18) * t19 + (-t123 * t74 + t2) * pkin(5)) * m(7) + (-t123 * t72 + t40) * pkin(5) + t23 * t99 - t18 * t96 - t22 * t97 + t19 * t98 + t416; -t275 * t96 + t123 * t98 + 0.2e1 * (t24 / 0.2e1 + t19 * t398 + t17 * t346) * m(7) + t15;];
tauc  = t1(:);

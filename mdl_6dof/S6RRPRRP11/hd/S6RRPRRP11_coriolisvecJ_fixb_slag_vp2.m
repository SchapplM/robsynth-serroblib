% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:18:27
% EndTime: 2018-11-23 17:18:39
% DurationCPUTime: 11.94s
% Computational Cost: add. (8130->597), mult. (18639->787), div. (0->0), fcn. (11496->6), ass. (0->274)
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
t366 = qJD(3) - t194;
t152 = qJD(2) * t245 + t366;
t92 = -t147 * t240 + t243 * t152;
t80 = -pkin(9) * t187 + t92;
t71 = pkin(4) * t220 + t80;
t93 = t147 * t243 + t152 * t240;
t81 = pkin(9) * t186 + t93;
t77 = t239 * t81;
t22 = t242 * t71 - t77;
t406 = qJ(6) * t123;
t18 = t22 - t406;
t293 = qJD(4) + qJD(5);
t211 = t233 + t293;
t17 = pkin(5) * t211 + t18;
t79 = t242 * t81;
t23 = t239 * t71 + t79;
t275 = t242 * t186 - t187 * t239;
t370 = qJ(6) * t275;
t19 = t23 + t370;
t295 = qJD(1) * qJD(2);
t279 = t244 * t295;
t216 = Ifges(7,3) * t279;
t217 = Ifges(6,3) * t279;
t340 = -t211 / 0.2e1;
t348 = t123 / 0.2e1;
t349 = -t123 / 0.2e1;
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
t356 = pkin(3) + pkin(7);
t197 = t356 * t302;
t177 = qJD(1) * t197;
t47 = -qJD(4) * t93 - t134 * t240 + t243 * t177;
t27 = pkin(4) * t279 - pkin(9) * t138 + t47;
t299 = qJD(4) * t243;
t300 = qJD(4) * t240;
t46 = t243 * t134 - t147 * t300 + t152 * t299 + t240 * t177;
t33 = pkin(9) * t139 + t46;
t6 = -qJD(5) * t23 - t239 * t33 + t242 * t27;
t2 = pkin(5) * t279 - qJ(6) * t56 - qJD(6) * t123 + t6;
t296 = qJD(5) * t242;
t297 = qJD(5) * t239;
t5 = t239 * t27 + t242 * t33 + t71 * t296 - t297 * t81;
t57 = -qJD(5) * t123 - t138 * t239 + t139 * t242;
t3 = qJ(6) * t57 + qJD(6) * t275 + t5;
t365 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t393 = Ifges(6,5) + Ifges(7,5);
t395 = Ifges(6,1) + Ifges(7,1);
t394 = Ifges(6,4) + Ifges(7,4);
t414 = t394 * t275;
t385 = t123 * t395 + t393 * t211 + t414;
t391 = Ifges(6,6) + Ifges(7,6);
t392 = Ifges(6,2) + Ifges(7,2);
t411 = t123 * t394;
t386 = t211 * t391 + t275 * t392 + t411;
t399 = -t275 / 0.2e1;
t51 = Ifges(7,6) * t57;
t52 = Ifges(6,6) * t57;
t53 = Ifges(7,5) * t56;
t54 = Ifges(6,5) * t56;
t74 = -pkin(5) * t275 + qJD(6) + t130;
t418 = t216 + t217 + t51 + t52 + t53 + t54 + t365 + (-t123 * t391 + t275 * t393) * t340 + (t123 * t19 + t17 * t275) * mrSges(7,3) + (t123 * t23 + t22 * t275) * mrSges(6,3) - t130 * (mrSges(6,1) * t123 + mrSges(6,2) * t275) - t74 * (mrSges(7,1) * t123 + mrSges(7,2) * t275) + t386 * t348 + (-t123 * t392 + t385 + t414) * t399 + (t395 * t275 - t411) * t349;
t228 = pkin(2) * t233;
t156 = qJD(1) * t262 + t228;
t113 = -t156 * t240 + t243 * t195;
t255 = -pkin(9) * t240 * t241 + pkin(4) * t244;
t332 = pkin(9) - t245;
t417 = -qJD(1) * t255 + t332 * t300 - t113;
t114 = t243 * t156 + t240 * t195;
t201 = t332 * t243;
t284 = t243 * t233;
t416 = pkin(9) * t284 + qJD(4) * t201 + t114;
t328 = Ifges(5,4) * t240;
t264 = Ifges(5,2) * t243 + t328;
t415 = -t264 / 0.2e1;
t398 = -qJD(1) / 0.2e1;
t310 = t242 * t243;
t127 = -t239 * t300 - t240 * t297 + t293 * t310;
t311 = t239 * t240;
t148 = -t233 * t311 + t242 * t284;
t308 = t127 + t148;
t257 = t239 * t243 + t242 * t240;
t128 = t293 * t257;
t251 = t257 * t241;
t149 = qJD(1) * t251;
t413 = t128 + t149;
t318 = t187 * Ifges(5,4);
t110 = t186 * Ifges(5,2) + t220 * Ifges(5,6) + t318;
t179 = Ifges(5,4) * t186;
t111 = Ifges(5,1) * t187 + Ifges(5,5) * t220 + t179;
t260 = t240 * t92 - t243 * t93;
t269 = mrSges(5,1) * t243 - mrSges(5,2) * t240;
t335 = -t243 / 0.2e1;
t336 = -t240 / 0.2e1;
t412 = t260 * mrSges(5,3) + t110 * t335 + t111 * t336 + t171 * t269;
t410 = qJD(3) + t227;
t408 = -t394 * t57 / 0.2e1 - t395 * t56 / 0.2e1 - t393 * t279 / 0.2e1;
t407 = -mrSges(3,1) + mrSges(4,2);
t200 = t332 * t240;
t384 = t200 * t297 - t201 * t296 + t239 * t417 - t242 * t416;
t132 = -t242 * t200 - t239 * t201;
t383 = -qJD(5) * t132 + t239 * t416 + t242 * t417;
t285 = -pkin(4) * t243 - pkin(3);
t368 = pkin(4) * t299 - t233 * t285 + t410;
t202 = -pkin(2) * t244 + t278;
t172 = t202 * qJD(1);
t205 = -t229 - t238;
t320 = Ifges(5,6) * t243;
t324 = Ifges(5,5) * t240;
t263 = t320 + t324;
t327 = Ifges(5,4) * t243;
t266 = Ifges(5,1) * t240 + t327;
t337 = t220 / 0.2e1;
t343 = t187 / 0.2e1;
t396 = qJD(2) / 0.2e1;
t397 = -qJD(2) / 0.2e1;
t405 = -t205 * mrSges(4,1) + t172 * mrSges(4,2) - Ifges(4,5) * t396 - Ifges(3,6) * t397 + t186 * t415 - t263 * t337 - t266 * t343 + t412 + ((-Ifges(3,2) - Ifges(4,3)) * t244 + (-Ifges(3,4) - Ifges(4,6)) * t241) * t398;
t256 = -t310 + t311;
t403 = t22 * t413 - t23 * t308 + t256 * t6 - t257 * t5;
t402 = t17 * t413 - t19 * t308 + t2 * t256 - t257 * t3;
t339 = t211 / 0.2e1;
t351 = t275 / 0.2e1;
t390 = t279 * t391 + t392 * t57 + t394 * t56;
t388 = -qJ(6) * t308 - qJD(6) * t257 + t384;
t387 = -pkin(5) * t305 + qJ(6) * t413 + qJD(6) * t256 + t383;
t369 = pkin(5) * t308 + t368;
t209 = t356 * t241;
t191 = t243 * t209;
t277 = pkin(9) * t244 - t178;
t106 = pkin(4) * t241 + t240 * t277 + t191;
t190 = t240 * t209;
t126 = t243 * t178 + t190;
t309 = t243 * t244;
t112 = -pkin(9) * t309 + t126;
t60 = t239 * t106 + t242 * t112;
t367 = t240 * t46 + t243 * t47;
t363 = t244 * t293;
t362 = t47 * mrSges(5,1) - t46 * mrSges(5,2) + Ifges(5,5) * t138 + Ifges(5,6) * t139;
t199 = -qJD(2) * pkin(2) + t410;
t226 = Ifges(3,4) * t305;
t287 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t288 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t290 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t322 = Ifges(4,6) * t244;
t361 = t290 * t123 + t287 * t211 + t288 * t275 + t17 * mrSges(7,1) + t199 * mrSges(4,1) + t22 * mrSges(6,1) + t92 * mrSges(5,1) + t220 * Ifges(5,3) + t187 * Ifges(5,5) + t186 * Ifges(5,6) + Ifges(3,1) * t233 / 0.2e1 + Ifges(3,5) * t396 + t226 / 0.2e1 + Ifges(4,4) * t397 + (-t241 * Ifges(4,2) - t322) * t398 - t172 * mrSges(4,3) - t19 * mrSges(7,2) - t23 * mrSges(6,2) - t93 * mrSges(5,2) + t391 * t351 + t393 * t348 + (Ifges(7,3) + Ifges(6,3)) * t339;
t360 = t56 / 0.2e1;
t359 = t57 / 0.2e1;
t355 = pkin(1) * mrSges(3,1);
t354 = pkin(1) * mrSges(3,2);
t345 = -t186 / 0.2e1;
t344 = -t187 / 0.2e1;
t338 = -t220 / 0.2e1;
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
t222 = t240 * pkin(4) + qJ(3);
t129 = -mrSges(5,1) * t186 + mrSges(5,2) * t187;
t207 = -mrSges(4,1) * t305 - qJD(2) * mrSges(4,3);
t306 = -t207 + t129;
t210 = t356 * t244;
t292 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t291 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t289 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t286 = m(4) * pkin(7) + mrSges(4,1);
t165 = pkin(4) * t309 + t210;
t15 = -t57 * mrSges(7,1) + t56 * mrSges(7,2);
t29 = -t239 * t80 - t79;
t59 = t242 * t106 - t112 * t239;
t232 = pkin(2) * t304;
t143 = t232 + t249;
t276 = -t143 * t240 + t243 * t197;
t131 = t200 * t239 - t242 * t201;
t196 = t356 * t304;
t273 = m(4) * t199 + (mrSges(4,1) + mrSges(3,3)) * t233 + t407 * qJD(2);
t272 = m(4) * t205 - mrSges(3,3) * t305 + t207 + t313;
t268 = mrSges(5,1) * t240 + mrSges(5,2) * t243;
t267 = Ifges(5,1) * t243 - t328;
t265 = -Ifges(5,2) * t240 + t327;
t118 = mrSges(5,1) * t279 - mrSges(5,3) * t138;
t119 = -mrSges(5,2) * t279 + mrSges(5,3) * t139;
t259 = t243 * t118 + t240 * t119;
t141 = -mrSges(5,2) * t220 + mrSges(5,3) * t186;
t142 = mrSges(5,1) * t220 - mrSges(5,3) * t187;
t258 = t243 * t141 - t240 * t142;
t48 = t255 * qJD(2) + (t243 * t277 - t190) * qJD(4) + t276;
t67 = t243 * t143 - t178 * t300 + t240 * t197 + t209 * t299;
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
t85 = t138 * Ifges(5,1) + t139 * Ifges(5,4) + Ifges(5,5) * t279;
t84 = t138 * Ifges(5,4) + t139 * Ifges(5,2) + Ifges(5,6) * t279;
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
t21 = t30 - t406;
t20 = t29 - t370;
t16 = -mrSges(6,1) * t57 + mrSges(6,2) * t56;
t8 = qJ(6) * t83 + qJD(6) * t157 + t9;
t7 = pkin(5) * t302 - qJ(6) * t82 + qJD(6) * t158 + t10;
t1 = [t210 * t88 + t160 * t192 - t196 * t129 + t165 * t16 + t67 * t141 + t68 * t142 + t133 * t73 + t130 * (-mrSges(6,1) * t83 + mrSges(6,2) * t82) + t125 * t118 + t126 * t119 + t115 * t15 + m(7) * (t115 * t24 + t17 * t7 + t19 * t8 + t2 * t36 + t3 * t37 + t58 * t74) + m(6) * (t10 * t22 + t102 * t165 + t130 * t133 + t23 * t9 + t5 * t60 + t59 * t6) + t102 * (-mrSges(6,1) * t157 - mrSges(6,2) * t158) + t24 * (-mrSges(7,1) * t157 - mrSges(7,2) * t158) + (t157 * t3 + t158 * t2 - t17 * t82 + t19 * t83) * mrSges(7,3) + (t157 * t5 + t158 * t6 - t22 * t82 + t23 * t83) * mrSges(6,3) + t158 * t408 + m(5) * (t125 * t47 + t126 * t46 + t153 * t210 - t171 * t196 + t67 * t93 + t68 * t92) + (-t145 * mrSges(4,3) + t53 / 0.2e1 + t51 / 0.2e1 + t216 / 0.2e1 + t54 / 0.2e1 + t52 / 0.2e1 + t217 / 0.2e1 + t218 / 0.2e1 + t288 * t57 + t290 * t56 + t362 + t365) * t241 + t7 * t98 + t10 * t99 + t8 * t96 + t9 * t97 + t74 * (-mrSges(7,1) * t83 + mrSges(7,2) * t82) + t58 * t72 + t59 * t41 + t60 * t43 + t37 * t42 + t36 * t40 + (t291 * qJD(2) + (-t202 * mrSges(4,2) + t241 * t289 - 0.2e1 * t355) * qJD(1) + t272 * pkin(7) - t405) * t304 + t385 * t82 / 0.2e1 + t386 * t83 / 0.2e1 + t390 * t157 / 0.2e1 + (t391 * t83 + t393 * t82) * t339 + (t157 * t392 - t158 * t394) * t359 + (t392 * t83 + t394 * t82) * t351 + (t157 * t394 - t158 * t395) * t360 + (t394 * t83 + t395 * t82) * t348 + (t139 * t415 + t145 * mrSges(4,2) + t153 * t269 - t138 * t266 / 0.2e1 + t85 * t336 + t84 * t335 - t286 * t198 + (t240 * t47 - t243 * t46) * mrSges(5,3) + (t111 * t335 + t240 * t110 / 0.2e1 + (t321 - t323) * t337 + t265 * t345 + t267 * t344 - t171 * t268 + (t240 * t93 + t243 * t92) * mrSges(5,3)) * qJD(4) + (t292 * qJD(2) + t273 * pkin(7) + ((-t324 / 0.2e1 - t320 / 0.2e1 - t289) * t244 - t202 * mrSges(4,3) - 0.2e1 * t354 - t290 * t158 + t288 * t157 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(5,3) / 0.2e1 + t286 * pkin(7) + t287) * t241) * qJD(1) + t361) * qJD(2)) * t244 + m(4) * (t145 * t202 + t160 * t172); (mrSges(7,1) * t308 - mrSges(7,2) * t413) * t74 + (mrSges(6,1) * t308 - mrSges(6,2) * t413) * t130 + t222 * t16 - t194 * t129 - t198 * mrSges(4,3) + t150 * t15 - t114 * t141 - t113 * t142 + t132 * t43 + t131 * t41 + (t259 + (-m(5) * t260 + t258) * qJD(4) + m(5) * t367) * t245 + t368 * t73 + (-m(4) * t172 - t192) * (-qJ(3) * t305 + t228) + m(4) * (-qJ(3) * t198 - qJD(3) * t205) + t369 * t72 - t367 * mrSges(5,3) + t256 * t408 + (t153 * qJ(3) - t113 * t92 - t114 * t93 + t171 * t366) * m(5) + t103 * t40 + t104 * t42 + qJ(3) * t88 - t390 * t257 / 0.2e1 + t102 * (mrSges(6,1) * t257 - mrSges(6,2) * t256) + t24 * (mrSges(7,1) * t257 - mrSges(7,2) * t256) + (-t256 * t394 - t257 * t392) * t359 + (-t256 * t395 - t257 * t394) * t360 + (t148 * t392 + t149 * t394) * t399 + t385 * (-t149 / 0.2e1 - t128 / 0.2e1) + (-t127 * t391 - t128 * t393) * t339 + (-t127 * t392 - t128 * t394) * t351 + (-t127 * t394 - t128 * t395) * t348 + t139 * t265 / 0.2e1 + t138 * t267 / 0.2e1 + t153 * t268 + t243 * t85 / 0.2e1 + t402 * mrSges(7,3) + t403 * mrSges(6,3) + t383 * t99 + t384 * t97 + (t102 * t222 + t130 * t368 + t131 * t6 + t132 * t5 + t22 * t383 + t23 * t384) * m(6) + t386 * (-t148 / 0.2e1 - t127 / 0.2e1) + t387 * t98 + t388 * t96 + (t103 * t2 + t104 * t3 + t150 * t24 + t17 * t387 + t19 * t388 + t369 * t74) * m(7) + (t148 * t391 + t149 * t393) * t340 + (t148 * t394 + t149 * t395) * t349 + t84 * t336 + (t263 * t338 + t264 * t345 + t266 * t344 + t412) * qJD(4) + (((t355 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t241) * qJD(1) + (-qJ(3) * mrSges(4,1) + t291) * qJD(2) + (-t272 + t313) * pkin(7) + t405) * t241 + ((t323 / 0.2e1 - t321 / 0.2e1 - pkin(2) * mrSges(4,1) - t290 * t256 - t288 * t257 + t292) * qJD(2) + (t354 - t322 / 0.2e1) * qJD(1) - t226 / 0.2e1 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t233 + ((-m(4) * pkin(2) + t407) * qJD(2) - t273) * pkin(7) - t361) * t244) * qJD(1) + t306 * qJD(3); -(t41 + t40) * t256 + t331 * t257 + t258 * qJD(4) + (t192 + t258) * t233 + (t286 * t305 - t306 - t72 - t73) * qJD(2) - m(4) * (-qJD(2) * t205 - t172 * t233) + t259 + t308 * t330 - t413 * t329 + (-qJD(2) * t74 - t402) * m(7) + (-qJD(2) * t130 - t403) * m(6) + (-qJD(2) * t171 - t220 * t260 + t367) * m(5); t225 * t40 - t171 * (mrSges(5,1) * t187 + mrSges(5,2) * t186) - t92 * t141 + t93 * t142 + (t186 * t92 + t187 * t93) * mrSges(5,3) - m(6) * (t22 * t29 + t23 * t30) + (-t187 * t73 + t242 * t41 + t331 * t239 + (-t239 * t329 + t242 * t330) * qJD(5) + (-t130 * t187 - t22 * t297 + t23 * t296 + t239 * t5 + t242 * t6) * m(6)) * pkin(4) + (-t17 * t20 - t19 * t21 - t74 * t90 + t2 * t225 + (-t17 * t297 + t19 * t296 + t239 * t3) * pkin(4)) * m(7) + t362 + (-Ifges(5,2) * t187 + t111 + t179) * t345 - t20 * t98 - t29 * t99 - t90 * t72 - t21 * t96 - t30 * t97 + t110 * t343 + (Ifges(5,1) * t186 - t318) * t344 + (Ifges(5,5) * t186 - Ifges(5,6) * t187) * t338 + t218 + t418; (-t123 * t72 + t40) * pkin(5) + t19 * t98 + t23 * t99 - t18 * t96 - t22 * t97 + (-(-t17 + t18) * t19 + (-t123 * t74 + t2) * pkin(5)) * m(7) + t418; -t275 * t96 + t123 * t98 + 0.2e1 * (t24 / 0.2e1 + t19 * t399 + t17 * t348) * m(7) + t15;];
tauc  = t1(:);

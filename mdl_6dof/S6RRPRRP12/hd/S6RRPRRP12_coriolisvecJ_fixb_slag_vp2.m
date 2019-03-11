% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP12
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:25
% EndTime: 2019-03-09 12:50:50
% DurationCPUTime: 11.58s
% Computational Cost: add. (8044->635), mult. (18220->829), div. (0->0), fcn. (11031->6), ass. (0->293)
t390 = Ifges(6,1) + Ifges(7,1);
t389 = Ifges(6,4) - Ifges(7,5);
t388 = Ifges(7,4) + Ifges(6,5);
t228 = sin(qJ(5));
t229 = sin(qJ(4));
t231 = cos(qJ(4));
t300 = qJD(5) * t228;
t302 = qJD(4) * t231;
t342 = cos(qJ(5));
t279 = t342 * qJD(5);
t377 = qJD(4) * t342 + t279;
t115 = -t228 * t302 - t229 * t377 - t231 * t300;
t230 = sin(qJ(2));
t246 = -t228 * t231 - t229 * t342;
t240 = t230 * t246;
t138 = qJD(1) * t240;
t312 = t115 + t138;
t303 = qJD(4) * t229;
t116 = -t228 * t303 - t229 * t300 + t231 * t377;
t222 = t230 * qJD(1);
t286 = t342 * t231;
t272 = t230 * t286;
t314 = t228 * t229;
t137 = qJD(1) * t272 - t222 * t314;
t311 = t116 + t137;
t232 = cos(qJ(2));
t309 = qJD(1) * t232;
t218 = pkin(7) * t309;
t177 = pkin(3) * t309 + t218;
t227 = qJD(2) * qJ(3);
t157 = t227 + t177;
t233 = -pkin(2) - pkin(8);
t276 = -qJ(3) * t230 - pkin(1);
t164 = t232 * t233 + t276;
t136 = t164 * qJD(1);
t216 = pkin(7) * t222;
t176 = -pkin(3) * t222 - t216;
t375 = qJD(3) - t176;
t140 = qJD(2) * t233 + t375;
t85 = -t136 * t229 + t140 * t231;
t86 = t231 * t136 + t229 * t140;
t256 = t229 * t85 - t231 * t86;
t265 = mrSges(5,1) * t231 - mrSges(5,2) * t229;
t343 = -t231 / 0.2e1;
t344 = -t229 / 0.2e1;
t209 = t222 + qJD(4);
t285 = t231 * t309;
t308 = qJD(2) * t229;
t242 = t285 + t308;
t306 = qJD(2) * t231;
t169 = -t229 * t309 + t306;
t328 = Ifges(5,4) * t169;
t98 = -t242 * Ifges(5,2) + Ifges(5,6) * t209 + t328;
t165 = Ifges(5,4) * t242;
t99 = Ifges(5,1) * t169 + Ifges(5,5) * t209 - t165;
t409 = -t256 * mrSges(5,3) - t157 * t265 - t343 * t98 - t344 * t99;
t408 = qJD(3) + t216;
t407 = Ifges(6,6) / 0.2e1;
t406 = -Ifges(7,6) / 0.2e1;
t109 = t169 * t228 + t242 * t342;
t307 = qJD(2) * t230;
t284 = t229 * t307;
t298 = qJD(2) * qJD(4);
t301 = qJD(4) * t232;
t127 = -t229 * t298 + (-t231 * t301 + t284) * qJD(1);
t283 = t229 * t301;
t241 = t230 * t306 + t283;
t128 = qJD(1) * t241 - t231 * t298;
t51 = -qJD(5) * t109 + t127 * t342 + t128 * t228;
t366 = t51 / 0.2e1;
t236 = t169 * t342 - t228 * t242;
t52 = qJD(5) * t236 + t228 * t127 - t128 * t342;
t365 = -t52 / 0.2e1;
t360 = pkin(3) + pkin(7);
t405 = -mrSges(3,1) + mrSges(4,2);
t75 = -pkin(9) * t242 + t86;
t288 = t342 * t75;
t74 = -pkin(9) * t169 + t85;
t25 = t228 * t74 + t288;
t404 = pkin(4) * t300 - t25;
t170 = -t286 + t314;
t145 = t170 * t232;
t287 = -pkin(4) * t231 - pkin(3);
t378 = pkin(4) * t302 - t222 * t287 + t408;
t184 = -pkin(2) * t232 + t276;
t158 = t184 * qJD(1);
t188 = -t218 - t227;
t190 = -mrSges(4,1) * t309 - qJD(2) * mrSges(4,3);
t320 = Ifges(5,6) * t231;
t323 = Ifges(5,5) * t229;
t259 = t320 + t323;
t326 = Ifges(5,4) * t231;
t262 = Ifges(5,1) * t229 + t326;
t345 = t209 / 0.2e1;
t350 = t169 / 0.2e1;
t391 = qJD(2) / 0.2e1;
t392 = -qJD(2) / 0.2e1;
t393 = -qJD(1) / 0.2e1;
t403 = (m(4) * t188 + qJD(2) * mrSges(3,2) - mrSges(3,3) * t309 + t190) * pkin(7) + t188 * mrSges(4,1) + Ifges(3,6) * t392 + (Ifges(3,4) * t230 + Ifges(3,2) * t232) * t393 + Ifges(4,5) * t391 + (-Ifges(4,6) * t230 - Ifges(4,3) * t232) * qJD(1) / 0.2e1 - t158 * mrSges(4,2) + t262 * t350 + t259 * t345 + t409;
t118 = pkin(4) * t242 + t157;
t195 = qJD(5) + t209;
t318 = t228 * t75;
t66 = pkin(4) * t209 + t74;
t21 = t342 * t66 - t318;
t379 = qJD(6) - t21;
t19 = -pkin(5) * t195 + t379;
t105 = Ifges(6,4) * t109;
t321 = Ifges(7,5) * t109;
t384 = t195 * t388 + t236 * t390 - t105 + t321;
t40 = pkin(5) * t109 - qJ(6) * t236 + t118;
t402 = mrSges(6,2) * t118 + mrSges(7,2) * t19 - mrSges(6,3) * t21 - mrSges(7,3) * t40 + t384 / 0.2e1;
t22 = t228 * t66 + t288;
t299 = qJD(1) * qJD(2);
t277 = t232 * t299;
t278 = t230 * t299;
t208 = pkin(2) * t278;
t258 = pkin(8) * t230 - qJ(3) * t232;
t304 = qJD(3) * t230;
t238 = qJD(2) * t258 - t304;
t122 = qJD(1) * t238 + t208;
t305 = qJD(2) * t232;
t179 = t360 * t305;
t163 = qJD(1) * t179;
t39 = -qJD(4) * t86 - t122 * t229 + t163 * t231;
t24 = pkin(4) * t277 - pkin(9) * t127 + t39;
t38 = t122 * t231 - t136 * t303 + t140 * t302 + t163 * t229;
t28 = pkin(9) * t128 + t38;
t5 = t228 * t24 + t279 * t66 + t28 * t342 - t300 * t75;
t6 = -qJD(5) * t22 - t228 * t28 + t24 * t342;
t401 = t170 * t6 - t21 * t312 - t22 * t311 + t246 * t5;
t2 = qJ(6) * t277 + qJD(6) * t195 + t5;
t20 = qJ(6) * t195 + t22;
t3 = -pkin(5) * t277 - t6;
t400 = -t170 * t3 + t19 * t312 + t2 * t246 - t20 * t311;
t67 = pkin(5) * t236 + qJ(6) * t109;
t399 = -t51 * Ifges(7,5) / 0.2e1 + t277 * t406 + Ifges(7,3) * t365;
t398 = Ifges(6,4) * t366 + Ifges(6,2) * t365 + t277 * t407;
t357 = -t109 / 0.2e1;
t356 = t109 / 0.2e1;
t347 = t195 / 0.2e1;
t353 = t236 / 0.2e1;
t327 = Ifges(5,4) * t229;
t260 = Ifges(5,2) * t231 + t327;
t396 = t260 / 0.2e1;
t395 = Ifges(5,2) * t344 + t326 / 0.2e1;
t387 = Ifges(6,6) - Ifges(7,6);
t386 = t277 * t388 - t389 * t52 + t390 * t51;
t385 = pkin(5) * t311 - qJ(6) * t312 + qJD(6) * t170 + t378;
t89 = -mrSges(7,2) * t109 + mrSges(7,3) * t195;
t330 = mrSges(6,3) * t109;
t90 = -mrSges(6,2) * t195 - t330;
t334 = t89 + t90;
t329 = mrSges(6,3) * t236;
t91 = mrSges(6,1) * t195 - t329;
t92 = -mrSges(7,1) * t195 + mrSges(7,2) * t236;
t333 = t91 - t92;
t146 = t246 * t232;
t192 = t360 * t230;
t172 = t229 * t192;
t114 = t164 * t231 + t172;
t376 = t229 * t38 + t231 * t39;
t374 = qJD(4) + qJD(5);
t373 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t313 = t231 * t232;
t100 = -pkin(9) * t313 + t114;
t173 = t231 * t192;
t275 = pkin(9) * t232 - t164;
t95 = pkin(4) * t230 + t229 * t275 + t173;
t332 = t100 * t342 + t228 * t95;
t253 = -pkin(9) * t229 * t230 + pkin(4) * t232;
t221 = pkin(2) * t307;
t132 = t221 + t238;
t274 = -t132 * t229 + t179 * t231;
t41 = t253 * qJD(2) + (t231 * t275 - t172) * qJD(4) + t274;
t63 = t132 * t231 - t164 * t303 + t179 * t229 + t192 * t302;
t49 = pkin(9) * t241 + t63;
t10 = -qJD(5) * t332 - t228 * t49 + t342 * t41;
t370 = t39 * mrSges(5,1) - t38 * mrSges(5,2) + Ifges(5,5) * t127 + Ifges(5,6) * t128;
t181 = -qJD(2) * pkin(2) + t408;
t215 = Ifges(3,4) * t309;
t290 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t291 = t407 + t406;
t293 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t368 = (m(4) * t181 + (mrSges(4,1) + mrSges(3,3)) * t222 + t405 * qJD(2)) * pkin(7) - t291 * t109 + t290 * t195 + t293 * t236 + t181 * mrSges(4,1) + t20 * mrSges(7,3) + t21 * mrSges(6,1) + t85 * mrSges(5,1) + Ifges(3,1) * t222 / 0.2e1 + Ifges(3,5) * t391 + t215 / 0.2e1 + Ifges(4,4) * t392 + (-t230 * Ifges(4,2) - Ifges(4,6) * t232) * t393 + Ifges(6,6) * t357 + Ifges(7,6) * t356 - Ifges(5,6) * t242 / 0.2e1 + Ifges(5,3) * t209 + Ifges(5,5) * t169 - t158 * mrSges(4,3) - t19 * mrSges(7,1) - t22 * mrSges(6,2) - t86 * mrSges(5,2) + t388 * t353 + (Ifges(6,3) + Ifges(7,2)) * t347;
t104 = Ifges(7,5) * t236;
t57 = Ifges(7,6) * t195 + Ifges(7,3) * t109 + t104;
t325 = Ifges(6,4) * t236;
t60 = -Ifges(6,2) * t109 + Ifges(6,6) * t195 + t325;
t367 = -mrSges(6,1) * t118 - mrSges(7,1) * t40 + mrSges(7,2) * t20 + mrSges(6,3) * t22 + t60 / 0.2e1 - t57 / 0.2e1;
t364 = t52 / 0.2e1;
t359 = pkin(1) * mrSges(3,1);
t358 = pkin(1) * mrSges(3,2);
t354 = -t236 / 0.2e1;
t351 = -t169 / 0.2e1;
t348 = -t195 / 0.2e1;
t346 = -t209 / 0.2e1;
t341 = pkin(4) * t169;
t340 = pkin(4) * t228;
t337 = pkin(9) - t233;
t34 = mrSges(6,1) * t277 - mrSges(6,3) * t51;
t35 = -mrSges(7,1) * t277 + mrSges(7,2) * t51;
t336 = -t34 + t35;
t33 = -mrSges(7,2) * t52 + mrSges(7,3) * t277;
t36 = -mrSges(6,2) * t277 - mrSges(6,3) * t52;
t335 = t36 + t33;
t217 = pkin(2) * t222;
t144 = qJD(1) * t258 + t217;
t101 = -t144 * t229 + t177 * t231;
t83 = qJD(1) * t253 + t101;
t102 = t144 * t231 + t177 * t229;
t88 = pkin(9) * t222 * t231 + t102;
t32 = t228 * t83 + t342 * t88;
t331 = mrSges(5,3) * t169;
t322 = Ifges(5,5) * t231;
t317 = t229 * Ifges(5,6);
t212 = pkin(4) * t229 + qJ(3);
t117 = mrSges(5,1) * t242 + mrSges(5,2) * t169;
t310 = -t190 + t117;
t193 = t360 * t232;
t297 = t342 * pkin(4);
t295 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t294 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t292 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t289 = m(4) * pkin(7) + mrSges(4,1);
t152 = pkin(4) * t313 + t193;
t183 = t337 * t231;
t273 = pkin(4) * t279;
t178 = t360 * t307;
t268 = t337 * t303;
t264 = mrSges(5,1) * t229 + mrSges(5,2) * t231;
t263 = Ifges(5,1) * t231 - t327;
t106 = mrSges(5,1) * t277 - mrSges(5,3) * t127;
t107 = -mrSges(5,2) * t277 + mrSges(5,3) * t128;
t255 = t106 * t231 + t107 * t229;
t239 = t242 * mrSges(5,3);
t130 = -mrSges(5,2) * t209 - t239;
t131 = mrSges(5,1) * t209 - t331;
t254 = t130 * t231 - t131 * t229;
t31 = -t228 * t88 + t342 * t83;
t53 = -t100 * t228 + t342 * t95;
t182 = t337 * t229;
t120 = -t182 * t342 - t183 * t228;
t247 = t182 * t228 - t183 * t342;
t9 = -t100 * t300 + t228 * t41 + t279 * t95 + t342 * t49;
t245 = t229 * t260;
t244 = t231 * t260;
t243 = -qJ(3) * t305 - t304;
t226 = qJD(2) * qJD(3);
t141 = -qJD(1) * t178 + t226;
t94 = -pkin(4) * t128 + t141;
t205 = Ifges(6,3) * t277;
t207 = Ifges(7,2) * t277;
t44 = Ifges(7,6) * t52;
t45 = Ifges(6,6) * t52;
t46 = Ifges(6,5) * t51;
t47 = Ifges(7,4) * t51;
t237 = t205 + t207 + t44 - t45 + t46 + t47 + t373;
t121 = -pkin(4) * t283 + (-pkin(7) + t287) * t307;
t214 = -t297 - pkin(5);
t211 = qJ(6) + t340;
t206 = Ifges(5,3) * t277;
t196 = t273 + qJD(6);
t180 = pkin(7) * t278 - t226;
t174 = (mrSges(4,2) * t232 - mrSges(4,3) * t230) * qJD(1);
t168 = qJD(4) * t183;
t148 = t221 + t243;
t134 = qJD(1) * t243 + t208;
t113 = -t164 * t229 + t173;
t103 = -pkin(5) * t246 + qJ(6) * t170 + t212;
t82 = -mrSges(5,1) * t128 + mrSges(5,2) * t127;
t80 = -pkin(5) * t145 - qJ(6) * t146 + t152;
t79 = Ifges(5,1) * t127 + Ifges(5,4) * t128 + Ifges(5,5) * t277;
t78 = Ifges(5,4) * t127 + Ifges(5,2) * t128 + Ifges(5,6) * t277;
t77 = qJD(2) * t272 - t146 * t374 - t228 * t284;
t76 = -qJD(2) * t240 + t145 * t374;
t71 = qJD(5) * t120 - t168 * t228 - t268 * t342;
t70 = qJD(5) * t247 - t168 * t342 + t228 * t268;
t69 = mrSges(6,1) * t109 + mrSges(6,2) * t236;
t68 = mrSges(7,1) * t109 - mrSges(7,3) * t236;
t64 = -qJD(4) * t114 + t274;
t56 = t341 + t67;
t50 = -pkin(5) * t230 - t53;
t48 = qJ(6) * t230 + t332;
t30 = -pkin(5) * t309 - t31;
t29 = qJ(6) * t309 + t32;
t26 = t342 * t74 - t318;
t18 = -pkin(5) * t77 - qJ(6) * t76 - qJD(6) * t146 + t121;
t17 = mrSges(6,1) * t52 + mrSges(6,2) * t51;
t16 = mrSges(7,1) * t52 - mrSges(7,3) * t51;
t11 = pkin(5) * t52 - qJ(6) * t51 - qJD(6) * t236 + t94;
t8 = -pkin(5) * t305 - t10;
t7 = qJ(6) * t305 + qJD(6) * t230 + t9;
t1 = [t193 * t82 + t48 * t33 + t50 * t35 + t53 * t34 + m(5) * (t113 * t39 + t114 * t38 + t141 * t193 - t157 * t178 + t63 * t86 + t64 * t85) + m(4) * (t134 * t184 + t148 * t158) + (t145 * t2 + t146 * t3) * mrSges(7,2) + (Ifges(6,2) * t357 - Ifges(7,3) * t356 + t347 * t387 + t353 * t389 + t367) * t77 + (-t127 * t262 / 0.2e1 - t128 * t260 / 0.2e1 + t141 * t265 + t134 * mrSges(4,2) + t79 * t344 + t78 * t343 - t289 * t180 + (t229 * t39 - t231 * t38) * mrSges(5,3) + (t285 * t395 + t263 * t351 - t157 * t264 + (t317 - t322) * t345 + t229 * t98 / 0.2e1 + t99 * t343 + (t229 * t86 + t231 * t85) * mrSges(5,3)) * qJD(4) + (t303 * t395 + (-t184 * mrSges(4,3) - 0.2e1 * t358 + (-t320 - t323 / 0.2e1 - t292) * t232 + t293 * t146 + t291 * t145 + (0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - t244 / 0.2e1 + t289 * pkin(7) + t290) * t230) * qJD(1) + (-t317 / 0.2e1 + t295) * qJD(2) + t368) * qJD(2)) * t232 + (t46 / 0.2e1 - t45 / 0.2e1 + t205 / 0.2e1 + t206 / 0.2e1 + t47 / 0.2e1 + t207 / 0.2e1 + t44 / 0.2e1 - t134 * mrSges(4,3) - t291 * t52 + t293 * t51 + t370 + t373) * t230 + m(7) * (t11 * t80 + t18 * t40 + t19 * t8 + t2 * t48 + t20 * t7 + t3 * t50) - t178 * t117 + t148 * t174 + t152 * t17 + t94 * (-mrSges(6,1) * t145 + mrSges(6,2) * t146) + t11 * (-mrSges(7,1) * t145 - mrSges(7,3) * t146) + t63 * t130 + t64 * t131 + t121 * t69 + m(6) * (t10 * t21 + t118 * t121 + t152 * t94 + t22 * t9 + t332 * t5 + t53 * t6) + t332 * t36 + t113 * t106 + t114 * t107 + t7 * t89 + t9 * t90 + t10 * t91 + t8 * t92 + t80 * t16 + (Ifges(6,4) * t357 + Ifges(7,5) * t356 + t388 * t347 + t390 * t353 + t402) * t76 + (t145 * t389 + t146 * t390) * t366 + t386 * t146 / 0.2e1 + ((-t245 / 0.2e1 + t294) * qJD(2) + (-mrSges(4,2) * t184 + t230 * t292 - 0.2e1 * t359) * qJD(1) + t403) * t307 + t145 * t398 + t145 * t399 + (Ifges(7,5) * t146 - Ifges(7,3) * t145) * t364 + (Ifges(6,4) * t146 + Ifges(6,2) * t145) * t365 + (t145 * t5 - t146 * t6) * mrSges(6,3) + t18 * t68; t231 * t79 / 0.2e1 + t212 * t17 + t310 * qJD(3) + (mrSges(6,1) * t311 + mrSges(6,2) * t312) * t118 + (mrSges(7,1) * t311 - mrSges(7,3) * t312) * t40 + (t141 * qJ(3) - t101 * t85 - t102 * t86 + t157 * t375) * m(5) + (t255 + (-m(5) * t256 + t254) * qJD(4) + m(5) * t376) * t233 + (-m(4) * t158 - t174) * (-qJ(3) * t309 + t217) + t378 * t69 + t127 * t263 / 0.2e1 + t141 * t264 - t333 * t71 + t334 * t70 + t335 * t120 - t176 * t117 - t180 * mrSges(4,3) - t102 * t130 - t101 * t131 + (t57 - t60) * (t137 / 0.2e1 + t116 / 0.2e1) + t103 * t16 - t29 * t89 - t32 * t90 - t31 * t91 - t30 * t92 + qJ(3) * t82 + t400 * mrSges(7,2) + t401 * mrSges(6,3) + m(4) * (-qJ(3) * t180 - qJD(3) * t188) + (t115 * t390 - t116 * t389) * t353 + t385 * t68 - t386 * t170 / 0.2e1 + (t259 * t346 + t262 * t351 + t308 * t396 - t409) * qJD(4) + (t247 * t6 + t120 * t5 + t212 * t94 + (-t32 + t70) * t22 + (-t31 - t71) * t21 + t378 * t118) * m(6) - t336 * t247 + (t103 * t11 - t247 * t3 + t120 * t2 + t385 * t40 + (-t29 + t70) * t20 + (-t30 + t71) * t19) * m(7) + (-Ifges(7,5) * t170 - Ifges(7,3) * t246) * t364 + (-Ifges(6,4) * t170 + Ifges(6,2) * t246) * t365 + t11 * (-mrSges(7,1) * t246 + mrSges(7,3) * t170) + t94 * (-mrSges(6,1) * t246 - mrSges(6,2) * t170) + (-t170 * t390 + t246 * t389) * t366 + t246 * t398 + t246 * t399 + (Ifges(6,4) * t115 - Ifges(7,5) * t138 - Ifges(6,2) * t116 - Ifges(7,3) * t137) * t357 + (-Ifges(6,4) * t138 + Ifges(7,5) * t115 + Ifges(6,2) * t137 + Ifges(7,3) * t116) * t356 + t384 * (t138 / 0.2e1 + t115 / 0.2e1) + (t137 * t387 - t138 * t388) * t348 + (t137 * t389 - t138 * t390) * t354 + (t115 * t388 - t116 * t387) * t347 + (((t359 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t230) * qJD(1) + (t245 / 0.2e1 + pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,1) + t294) * qJD(2) - t403) * t230 + (-t215 / 0.2e1 + (t244 / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t222 + t302 * t396 + (t358 + (t320 / 0.2e1 - Ifges(4,6) / 0.2e1) * t232) * qJD(1) + (t322 / 0.2e1 - pkin(2) * mrSges(4,1) + t291 * t246 - t293 * t170 + (-m(4) * pkin(2) + t405) * pkin(7) + t295) * qJD(2) - t368) * t232) * qJD(1) - t376 * mrSges(5,3) + t128 * t395 + t78 * t344; -t335 * t246 + t336 * t170 + t254 * qJD(4) + (t174 + t254) * t222 + (t289 * t309 - t310 - t68 - t69) * qJD(2) - m(4) * (-qJD(2) * t188 - t158 * t222) + t255 + t311 * t334 + t312 * t333 + (-qJD(2) * t40 - t400) * m(7) + (-qJD(2) * t118 - t401) * m(6) + (-qJD(2) * t157 - t209 * t256 + t376) * m(5); t211 * t33 + t214 * t35 + t196 * t89 + t206 + (-t130 - t239) * t85 - t157 * (t169 * mrSges(5,1) - mrSges(5,2) * t242) + t237 + t90 * t273 - t69 * t341 + (t2 * t211 + t214 * t3 - t40 * t56 + (t196 - t26) * t20 + t404 * t19) * m(7) + (-t118 * t341 + t21 * t25 - t22 * t26 + (t342 * t6 + t228 * t5 + (-t21 * t228 + t22 * t342) * qJD(5)) * pkin(4)) * m(6) + (t131 + t331) * t86 - t334 * t26 + t370 + (-Ifges(6,2) * t356 + Ifges(7,3) * t357 - t387 * t348 - t389 * t354 + t367) * t236 - t56 * t68 + t34 * t297 + t36 * t340 + (-Ifges(5,5) * t242 - Ifges(5,6) * t169) * t346 + t98 * t350 + (-Ifges(5,1) * t242 - t328) * t351 + (-Ifges(5,2) * t169 - t165 + t99) * t242 / 0.2e1 - t404 * t333 + (-Ifges(6,4) * t356 - Ifges(7,5) * t357 - t388 * t348 - t390 * t354 + t402) * t109; t237 - pkin(5) * t35 + (t329 + t333) * t22 + (-t330 - t334) * t21 - t118 * (mrSges(6,1) * t236 - mrSges(6,2) * t109) + t60 * t353 - t40 * (mrSges(7,1) * t236 + mrSges(7,3) * t109) + (Ifges(7,3) * t236 - t321) * t357 + qJD(6) * t89 + qJ(6) * t33 + (t109 * t19 + t20 * t236) * mrSges(7,2) - t67 * t68 + (-t109 * t388 - t236 * t387) * t348 + (-pkin(5) * t3 + qJ(6) * t2 - t19 * t22 + t20 * t379 - t40 * t67) * m(7) + (-Ifges(6,2) * t236 - t105 + t384) * t356 + (-t109 * t390 + t104 - t325 + t57) * t354; t236 * t68 - t195 * t89 + 0.2e1 * (t3 / 0.2e1 + t40 * t353 + t20 * t348) * m(7) + t35;];
tauc  = t1(:);

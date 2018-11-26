% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2018-11-23 16:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:56:38
% EndTime: 2018-11-23 16:56:52
% DurationCPUTime: 14.03s
% Computational Cost: add. (9525->595), mult. (24987->766), div. (0->0), fcn. (18117->8), ass. (0->287)
t233 = sin(qJ(2));
t235 = cos(qJ(2));
t303 = sin(pkin(9));
t305 = cos(pkin(9));
t211 = t305 * t233 + t303 * t235;
t197 = t211 * qJD(1);
t232 = sin(qJ(4));
t234 = cos(qJ(4));
t168 = qJD(2) * t234 - t197 * t232;
t169 = qJD(2) * t232 + t197 * t234;
t231 = sin(pkin(10));
t304 = cos(pkin(10));
t106 = -t304 * t168 + t169 * t231;
t357 = -t106 / 0.2e1;
t237 = -t233 * t303 + t235 * t305;
t196 = t237 * qJD(1);
t191 = qJD(4) - t196;
t343 = -t191 / 0.2e1;
t342 = t191 / 0.2e1;
t243 = t231 * t168 + t169 * t304;
t354 = -t243 / 0.2e1;
t353 = t243 / 0.2e1;
t401 = Ifges(7,2) + Ifges(6,3);
t425 = Ifges(5,3) + t401;
t286 = qJD(4) * t232;
t297 = t196 * t232;
t424 = t286 - t297;
t324 = -qJ(3) - pkin(7);
t220 = t324 * t235;
t214 = qJD(1) * t220;
t201 = t303 * t214;
t219 = t324 * t233;
t213 = qJD(1) * t219;
t206 = qJD(2) * pkin(2) + t213;
t157 = t206 * t305 + t201;
t150 = -qJD(2) * pkin(3) - t157;
t104 = -t168 * pkin(4) + qJD(5) + t150;
t229 = -pkin(2) * t235 - pkin(1);
t291 = qJD(1) * t229;
t215 = qJD(3) + t291;
t129 = -t196 * pkin(3) - t197 * pkin(8) + t215;
t267 = t305 * t214;
t158 = t303 * t206 - t267;
t151 = qJD(2) * pkin(8) + t158;
t84 = t234 * t129 - t151 * t232;
t67 = -qJ(5) * t169 + t84;
t53 = pkin(4) * t191 + t67;
t85 = t129 * t232 + t151 * t234;
t68 = qJ(5) * t168 + t85;
t62 = t304 * t68;
t18 = t231 * t53 + t62;
t14 = qJ(6) * t191 + t18;
t41 = t106 * pkin(5) - qJ(6) * t243 + t104;
t423 = Ifges(6,4) * t353 + Ifges(7,5) * t354 + Ifges(6,6) * t342 + Ifges(7,6) * t343 + (Ifges(6,2) + Ifges(7,3)) * t357 - t104 * mrSges(6,1) - t41 * mrSges(7,1) + t14 * mrSges(7,2) + t18 * mrSges(6,3);
t402 = Ifges(7,4) + Ifges(6,5);
t404 = Ifges(6,1) + Ifges(7,1);
t422 = t402 * t342 + t404 * t353;
t200 = t237 * qJD(2);
t186 = qJD(1) * t200;
t123 = qJD(4) * t168 + t186 * t234;
t124 = -qJD(4) * t169 - t186 * t232;
t76 = t123 * t231 - t124 * t304;
t365 = -t76 / 0.2e1;
t77 = t123 * t304 + t231 * t124;
t363 = t77 / 0.2e1;
t198 = t211 * qJD(2);
t185 = qJD(1) * t198;
t344 = t185 / 0.2e1;
t403 = Ifges(6,4) - Ifges(7,5);
t400 = Ifges(6,6) - Ifges(7,6);
t341 = -t196 / 0.2e1;
t421 = Ifges(4,2) * t341;
t210 = t231 * t234 + t304 * t232;
t130 = t210 * t196;
t195 = t210 * qJD(4);
t391 = t130 - t195;
t264 = t304 * t234;
t242 = -t231 * t232 + t264;
t131 = t242 * t196;
t199 = t242 * qJD(4);
t420 = -t199 + t131;
t356 = t106 / 0.2e1;
t418 = Ifges(6,2) * t357 - Ifges(7,3) * t356 + t400 * t342 + t403 * t353 + t423;
t308 = t231 * t68;
t17 = t304 * t53 - t308;
t13 = -t191 * pkin(5) + qJD(6) - t17;
t417 = t104 * mrSges(6,2) + mrSges(7,2) * t13 - mrSges(6,3) * t17 - t41 * mrSges(7,3) + t357 * t403 + t422;
t416 = mrSges(4,2) * t197;
t159 = t213 * t303 - t267;
t389 = pkin(4) * t424 - t159;
t268 = qJD(2) * t324;
t193 = qJD(3) * t235 + t233 * t268;
t176 = t193 * qJD(1);
t194 = -t233 * qJD(3) + t235 * t268;
t236 = qJD(1) * t194;
t122 = t176 * t305 + t236 * t303;
t284 = qJD(1) * qJD(2);
t271 = t233 * t284;
t258 = pkin(2) * t271;
t126 = pkin(3) * t185 - pkin(8) * t186 + t258;
t38 = -qJD(4) * t85 - t122 * t232 + t234 * t126;
t12 = pkin(4) * t185 - qJ(5) * t123 - qJD(5) * t169 + t38;
t285 = qJD(4) * t234;
t37 = t234 * t122 + t232 * t126 + t129 * t285 - t151 * t286;
t16 = qJ(5) * t124 + qJD(5) * t168 + t37;
t3 = t12 * t304 - t231 * t16;
t2 = -t185 * pkin(5) - t3;
t364 = t76 / 0.2e1;
t121 = t176 * t303 - t305 * t236;
t83 = -t124 * pkin(4) + t121;
t9 = t76 * pkin(5) - t77 * qJ(6) - qJD(6) * t243 + t83;
t413 = mrSges(6,2) * t83 + mrSges(7,2) * t2 - mrSges(6,3) * t3 - mrSges(7,3) * t9 + Ifges(7,5) * t364 + 0.2e1 * t344 * t402 + 0.2e1 * t363 * t404 + (Ifges(6,4) + t403) * t365;
t412 = -Ifges(6,2) * t356 + Ifges(7,3) * t357 - t354 * t403 + t423;
t410 = -Ifges(6,4) * t356 - Ifges(7,5) * t357 - t402 * t343 - t354 * t404 + t417;
t408 = Ifges(6,4) * t357 + Ifges(7,5) * t356 + t417 + t422;
t407 = -Ifges(4,6) / 0.2e1;
t406 = -Ifges(6,6) / 0.2e1;
t288 = qJD(2) * t233;
t405 = pkin(2) * t288;
t54 = -mrSges(7,2) * t76 + mrSges(7,3) * t185;
t55 = -mrSges(6,2) * t185 - mrSges(6,3) * t76;
t397 = t54 + t55;
t56 = mrSges(6,1) * t185 - mrSges(6,3) * t77;
t57 = -t185 * mrSges(7,1) + t77 * mrSges(7,2);
t396 = t57 - t56;
t93 = -mrSges(7,2) * t106 + mrSges(7,3) * t191;
t94 = -mrSges(6,2) * t191 - mrSges(6,3) * t106;
t323 = t93 + t94;
t95 = mrSges(6,1) * t191 - mrSges(6,3) * t243;
t96 = -mrSges(7,1) * t191 + mrSges(7,2) * t243;
t322 = t96 - t95;
t395 = t215 * mrSges(4,2);
t394 = -pkin(5) * t391 + qJ(6) * t420 - qJD(6) * t210 + t389;
t393 = Ifges(4,5) * qJD(2);
t311 = t197 * mrSges(4,3);
t292 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t168 + mrSges(5,2) * t169 + t311;
t392 = Ifges(5,5) * t123 + Ifges(5,6) * t124;
t156 = -pkin(3) * t237 - t211 * pkin(8) + t229;
t165 = t219 * t303 - t220 * t305;
t161 = t234 * t165;
t103 = t232 * t156 + t161;
t388 = -t232 * t38 + t234 * t37;
t289 = qJD(1) * t235;
t301 = Ifges(3,6) * qJD(2);
t318 = Ifges(3,4) * t233;
t387 = pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t289) + t301 / 0.2e1 + (t235 * Ifges(3,2) + t318) * qJD(1) / 0.2e1;
t386 = t185 * t425 - t400 * t76 + t402 * t77 + t392;
t385 = qJD(2) * t407;
t309 = t197 * Ifges(4,4);
t384 = -t309 / 0.2e1 + t421;
t256 = mrSges(5,1) * t232 + mrSges(5,2) * t234;
t244 = t150 * t256;
t251 = Ifges(5,5) * t234 - Ifges(5,6) * t232;
t316 = Ifges(5,4) * t234;
t253 = -Ifges(5,2) * t232 + t316;
t317 = Ifges(5,4) * t232;
t255 = Ifges(5,1) * t234 - t317;
t335 = t234 / 0.2e1;
t336 = -t232 / 0.2e1;
t345 = t169 / 0.2e1;
t314 = t169 * Ifges(5,4);
t98 = t168 * Ifges(5,2) + t191 * Ifges(5,6) + t314;
t166 = Ifges(5,4) * t168;
t99 = t169 * Ifges(5,1) + t191 * Ifges(5,5) + t166;
t381 = t251 * t342 + t255 * t345 + t168 * t253 / 0.2e1 + t244 + t99 * t335 + t98 * t336;
t4 = t231 * t12 + t304 * t16;
t1 = qJ(6) * t185 + qJD(6) * t191 + t4;
t380 = t38 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t37 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t375 = t215 * mrSges(4,1) + t84 * mrSges(5,1) + t17 * mrSges(6,1) - t13 * mrSges(7,1) - t85 * mrSges(5,2) - t18 * mrSges(6,2) + t14 * mrSges(7,3) + t384 + t385;
t373 = mrSges(6,1) * t83 + mrSges(7,1) * t9 - mrSges(7,2) * t1 - mrSges(6,3) * t4 + 0.2e1 * Ifges(7,3) * t364 - t77 * Ifges(6,4) / 0.2e1 + t185 * t406 + (-t403 + Ifges(7,5)) * t363 + (-t400 + Ifges(7,6)) * t344 + (-t365 + t364) * Ifges(6,2);
t372 = -0.2e1 * pkin(1);
t370 = Ifges(5,3) / 0.2e1;
t352 = t123 / 0.2e1;
t351 = t124 / 0.2e1;
t347 = -t168 / 0.2e1;
t346 = -t169 / 0.2e1;
t332 = pkin(4) * t169;
t331 = pkin(4) * t231;
t290 = qJD(1) * t233;
t330 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t290);
t245 = -qJ(5) * t200 - qJD(5) * t211;
t139 = t193 * t305 + t194 * t303;
t141 = pkin(3) * t198 - pkin(8) * t200 + t405;
t260 = -t139 * t232 + t234 * t141;
t300 = qJ(5) * t211;
t28 = pkin(4) * t198 + t245 * t234 + (-t161 + (-t156 + t300) * t232) * qJD(4) + t260;
t272 = t211 * t285;
t276 = t234 * t139 + t232 * t141 + t156 * t285;
t32 = -qJ(5) * t272 + (-qJD(4) * t165 + t245) * t232 + t276;
t8 = t231 * t28 + t304 * t32;
t296 = t196 * t234;
t281 = pkin(2) * t290;
t140 = pkin(3) * t197 - pkin(8) * t196 + t281;
t160 = t213 * t305 + t201;
t90 = t234 * t140 - t160 * t232;
t64 = pkin(4) * t197 - qJ(5) * t296 + t90;
t91 = t232 * t140 + t234 * t160;
t82 = -qJ(5) * t297 + t91;
t30 = t231 * t64 + t304 * t82;
t102 = t234 * t156 - t165 * t232;
t81 = -pkin(4) * t237 - t234 * t300 + t102;
t295 = t211 * t232;
t87 = -qJ(5) * t295 + t103;
t40 = t231 * t81 + t304 * t87;
t321 = mrSges(4,3) * t196;
t320 = mrSges(5,3) * t168;
t319 = mrSges(5,3) * t169;
t190 = Ifges(4,4) * t196;
t315 = t168 * Ifges(5,6);
t313 = t169 * Ifges(5,5);
t310 = t197 * Ifges(4,1);
t302 = Ifges(3,5) * qJD(2);
t164 = -t305 * t219 - t220 * t303;
t299 = t121 * t164;
t273 = t303 * pkin(2);
t226 = t273 + pkin(8);
t293 = qJ(5) + t226;
t287 = qJD(4) * t211;
t275 = t305 * pkin(2);
t274 = t304 * pkin(4);
t34 = t76 * mrSges(6,1) + t77 * mrSges(6,2);
t33 = t76 * mrSges(7,1) - t77 * mrSges(7,3);
t270 = t235 * t284;
t262 = t185 * mrSges(4,1) + t186 * mrSges(4,2);
t261 = t293 * t232;
t259 = qJD(4) * t293;
t228 = -t275 - pkin(3);
t257 = mrSges(5,1) * t234 - mrSges(5,2) * t232;
t254 = Ifges(5,1) * t232 + t316;
t252 = Ifges(5,2) * t234 + t317;
t250 = Ifges(5,5) * t232 + Ifges(5,6) * t234;
t249 = -t232 * t37 - t234 * t38;
t248 = -t232 * t85 - t234 * t84;
t247 = t232 * t84 - t234 * t85;
t138 = t193 * t303 - t305 * t194;
t127 = -mrSges(5,2) * t191 + t320;
t128 = mrSges(5,1) * t191 - t319;
t246 = t127 * t234 - t128 * t232;
t137 = pkin(4) * t295 + t164;
t7 = -t231 * t32 + t28 * t304;
t29 = -t231 * t82 + t304 * t64;
t39 = -t231 * t87 + t304 * t81;
t216 = -t234 * pkin(4) + t228;
t92 = t138 + (t200 * t232 + t272) * pkin(4);
t238 = -qJD(5) * t232 - t234 * t259;
t230 = Ifges(3,4) * t289;
t227 = -t274 - pkin(5);
t223 = qJ(6) + t331;
t207 = t293 * t234;
t205 = Ifges(3,1) * t290 + t230 + t302;
t174 = -qJD(2) * mrSges(4,2) + t321;
t173 = qJD(5) * t234 - t232 * t259;
t154 = -mrSges(4,1) * t196 + t416;
t153 = t207 * t304 - t231 * t261;
t152 = t207 * t231 + t261 * t304;
t147 = t190 + t310 + t393;
t145 = t242 * t211;
t144 = t210 * t211;
t142 = -pkin(5) * t242 - t210 * qJ(6) + t216;
t115 = t173 * t304 + t231 * t238;
t114 = t173 * t231 - t238 * t304;
t101 = -mrSges(5,2) * t185 + mrSges(5,3) * t124;
t100 = mrSges(5,1) * t185 - mrSges(5,3) * t123;
t97 = t191 * Ifges(5,3) + t313 + t315;
t89 = t200 * t242 - t210 * t287;
t88 = t211 * t231 * t286 - t200 * t210 - t264 * t287;
t80 = -mrSges(5,1) * t124 + mrSges(5,2) * t123;
t66 = mrSges(6,1) * t106 + mrSges(6,2) * t243;
t65 = mrSges(7,1) * t106 - mrSges(7,3) * t243;
t60 = t123 * Ifges(5,1) + t124 * Ifges(5,4) + t185 * Ifges(5,5);
t59 = t123 * Ifges(5,4) + t124 * Ifges(5,2) + t185 * Ifges(5,6);
t58 = t144 * pkin(5) - t145 * qJ(6) + t137;
t49 = Ifges(7,4) * t243 + t191 * Ifges(7,2) + t106 * Ifges(7,6);
t48 = Ifges(6,5) * t243 - t106 * Ifges(6,6) + t191 * Ifges(6,3);
t45 = pkin(5) * t243 + qJ(6) * t106 + t332;
t43 = -qJD(4) * t103 + t260;
t42 = -t165 * t286 + t276;
t36 = pkin(5) * t237 - t39;
t35 = -qJ(6) * t237 + t40;
t22 = t304 * t67 - t308;
t21 = t231 * t67 + t62;
t20 = -t197 * pkin(5) - t29;
t19 = qJ(6) * t197 + t30;
t11 = -t88 * pkin(5) - t89 * qJ(6) - t145 * qJD(6) + t92;
t6 = -t198 * pkin(5) - t7;
t5 = qJ(6) * t198 - qJD(6) * t237 + t8;
t10 = [t418 * t88 + (-t157 * t200 - t158 * t198 - t165 * t185) * mrSges(4,3) + t373 * t144 + m(7) * (t1 * t35 + t11 * t41 + t13 * t6 + t14 * t5 + t2 * t36 + t58 * t9) + m(6) * (t104 * t92 + t137 * t83 + t17 * t7 + t18 * t8 + t3 * t39 + t4 * t40) + t408 * t89 + (-(t370 + Ifges(4,2)) * t185 - t386 / 0.2e1 - t392 / 0.2e1 - mrSges(4,1) * qJD(1) * t405 + mrSges(4,3) * t122 + Ifges(4,4) * t186 - Ifges(6,6) * t365 - Ifges(7,6) * t364 - t401 * t344 - t402 * t363 - t380) * t237 + (mrSges(4,3) * t186 + t80) * t164 + t413 * t145 + (t255 * t352 + t253 * t351 + t251 * t344 + t59 * t336 + t60 * t335 + Ifges(4,1) * t186 - Ifges(4,4) * t185 + (mrSges(4,3) + t256) * t121 + t249 * mrSges(5,3) + (t150 * t257 + t252 * t347 + t254 * t346 + t250 * t343 + t99 * t336 - t234 * t98 / 0.2e1 + t247 * mrSges(5,3)) * qJD(4)) * t211 + (-t301 / 0.2e1 + (mrSges(3,1) * t372 - 0.3e1 / 0.2e1 * t318 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t235) * qJD(1) + (t416 + t154 + m(4) * (t215 + t291)) * pkin(2) - t387) * t288 + (t147 / 0.2e1 + t190 / 0.2e1 + t310 / 0.2e1 + t395 + t248 * mrSges(5,3) + t381) * t200 + ((Ifges(7,6) / 0.2e1 + t406) * t106 + t313 / 0.2e1 + t315 / 0.2e1 + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t243 + t97 / 0.2e1 + t49 / 0.2e1 + t48 / 0.2e1 + t375 + (t370 + Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t191 + t384) * t198 + (Ifges(4,5) * t200 / 0.2e1 + t198 * t407 + (-t330 + t205 / 0.2e1 + t302 / 0.2e1 + (mrSges(3,2) * t372 + 0.3e1 / 0.2e1 * Ifges(3,4) * t235) * qJD(1)) * t235) * qJD(2) + t229 * t262 + t35 * t54 + t40 * t55 + t39 * t56 + t36 * t57 + t58 * t33 + t11 * t65 + t92 * t66 + t5 * t93 + t8 * t94 + t7 * t95 + t6 * t96 + t102 * t100 + t103 * t101 + t42 * t127 + t43 * t128 + t292 * t138 + t137 * t34 + m(5) * (t102 * t38 + t103 * t37 + t138 * t150 + t42 * t85 + t43 * t84 + t299) + m(4) * (t122 * t165 - t138 * t157 + t139 * t158 + t299) + t139 * t174; -m(5) * (t84 * t90 + t85 * t91) - t418 * t195 + (t190 + t147) * t341 + ((-t121 * t305 + t122 * t303) * pkin(2) + t157 * t159 - t158 * t160 - t215 * t281) * m(4) - (-Ifges(3,2) * t290 + t205 + t230) * t289 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t233 + mrSges(3,2) * t235) - t233 * (Ifges(3,1) * t235 - t318) / 0.2e1) * qJD(1) ^ 2 + (m(5) * t228 - mrSges(4,1) - t257) * t121 + t387 * t290 + (m(5) * t388 - t232 * t100 + t234 * t101 - t127 * t286 - t128 * t285) * t226 + (-t152 * t3 + t153 * t4 + t216 * t83 + (-t30 + t115) * t18 + (-t29 - t114) * t17 + t389 * t104) * m(6) - t410 * t131 + (-t400 * t130 + t196 * t251 + t197 * t425) * t343 + t389 * t66 + (-m(5) * t150 - t292) * t159 + Ifges(3,5) * t270 + t408 * t199 + (-mrSges(3,1) * t270 + mrSges(3,2) * t271) * pkin(7) - t373 * t242 + (m(5) * t226 * t248 + t381) * qJD(4) - (-Ifges(5,5) * t346 - Ifges(5,6) * t347 - Ifges(6,6) * t356 - Ifges(7,6) * t357 - t354 * t402 + t375 + t385 + t421) * t197 + t412 * t130 + t413 * t210 + (-t424 * t85 + (-t285 + t296) * t84 + t388) * mrSges(5,3) + (t1 * t153 + t142 * t9 + t152 * t2 + t394 * t41 + (-t19 + t115) * t14 + (-t20 + t114) * t13) * m(7) + t394 * t65 + (-t244 + t255 * t346 + t253 * t347 - t393 / 0.2e1 - t395) * t196 + (-t185 * t273 - t186 * t275) * mrSges(4,3) - (Ifges(4,1) * t196 - t309 + t48 + t49 + t97) * t197 / 0.2e1 + t322 * t114 + t323 * t115 + t396 * t152 + t397 * t153 + t157 * t321 + t289 * t330 + t59 * t335 + t250 * t344 + t252 * t351 + t254 * t352 - Ifges(3,6) * t271 - t19 * t93 - t30 * t94 - t29 * t95 - t20 * t96 - t154 * t281 - (Ifges(3,5) * t235 - Ifges(3,6) * t233) * t284 / 0.2e1 - t122 * mrSges(4,2) - t91 * t127 - t90 * t128 - t99 * t296 / 0.2e1 + t142 * t33 + t98 * t297 / 0.2e1 - t160 * t174 + t158 * t311 - Ifges(4,6) * t185 + Ifges(4,5) * t186 + t216 * t34 + t228 * t80 + t232 * t60 / 0.2e1; t234 * t100 + t232 * t101 + t397 * t210 - t396 * t242 + t246 * qJD(4) + (-t174 - t246) * t196 - (t65 + t66 + t292) * t197 + t262 - t420 * t323 - t391 * t322 + (t1 * t210 - t13 * t391 - t14 * t420 - t197 * t41 - t2 * t242) * m(7) + (-t104 * t197 + t17 * t391 - t18 * t420 + t210 * t4 + t242 * t3) * m(6) + (-t150 * t197 - t191 * t247 - t249) * m(5) + (t157 * t197 - t158 * t196 + t258) * m(4); -t323 * t22 - t322 * t21 + t380 + t386 + ((t231 * t4 + t3 * t304) * pkin(4) - t104 * t332 + t17 * t21 - t18 * t22) * m(6) + (-t343 * t400 + t412) * t243 + (t1 * t223 - t13 * t21 + t2 * t227 - t41 * t45 + (qJD(6) - t22) * t14) * m(7) + (-Ifges(5,2) * t169 + t166 + t99) * t347 + (t320 - t127) * t84 + (t319 + t128) * t85 + t56 * t274 + t55 * t331 + (Ifges(5,5) * t168 - Ifges(5,6) * t169) * t343 + t98 * t345 + (Ifges(5,1) * t168 - t314) * t346 - t45 * t65 + qJD(6) * t93 - t150 * (mrSges(5,1) * t169 + mrSges(5,2) * t168) + t223 * t54 + t227 * t57 - t66 * t332 + t410 * t106; -t322 * t243 + t323 * t106 + t33 + t34 + (t106 * t14 - t13 * t243 + t9) * m(7) + (t106 * t18 + t17 * t243 + t83) * m(6); t243 * t65 - t191 * t93 + 0.2e1 * (t2 / 0.2e1 + t41 * t353 + t14 * t343) * m(7) + t57;];
tauc  = t10(:);

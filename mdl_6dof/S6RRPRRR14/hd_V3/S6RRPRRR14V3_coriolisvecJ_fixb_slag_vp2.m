% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR14V3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:37
% EndTime: 2019-04-12 15:04:06
% DurationCPUTime: 10.55s
% Computational Cost: add. (5041->604), mult. (13325->828), div. (0->0), fcn. (9489->8), ass. (0->280)
t123 = sin(qJ(4));
t124 = sin(qJ(2));
t127 = cos(qJ(4));
t215 = qJD(1) * t127;
t105 = qJD(2) * t123 + t124 * t215;
t128 = cos(qJ(2));
t214 = qJD(1) * t128;
t116 = qJD(4) - t214;
t122 = sin(qJ(5));
t126 = cos(qJ(5));
t72 = -t105 * t122 + t116 * t126;
t69 = qJD(6) - t72;
t278 = Ifges(7,3) * t69;
t216 = qJD(1) * t124;
t195 = t123 * t216;
t212 = qJD(2) * t127;
t104 = -t195 + t212;
t103 = qJD(5) - t104;
t121 = sin(qJ(6));
t125 = cos(qJ(6));
t73 = t105 * t126 + t116 * t122;
t45 = t103 * t125 - t121 * t73;
t279 = Ifges(7,6) * t45;
t46 = t103 * t121 + t125 * t73;
t281 = Ifges(7,5) * t46;
t17 = t278 + t279 + t281;
t249 = Ifges(6,6) * t103;
t280 = Ifges(6,2) * t72;
t283 = Ifges(6,4) * t73;
t30 = t249 + t280 + t283;
t369 = t30 / 0.2e1 - t17 / 0.2e1;
t182 = qJ(3) * t195;
t101 = qJ(3) * t212 - t182;
t82 = qJD(3) * t122 + t101 * t126;
t382 = -t82 * mrSges(6,3) - t369;
t381 = Ifges(5,2) / 0.2e1;
t97 = t105 * qJ(3);
t47 = -t121 * t82 + t125 * t97;
t48 = t121 * t97 + t125 * t82;
t328 = t121 * t48 + t125 * t47;
t380 = mrSges(7,3) * t328;
t213 = qJD(2) * t124;
t188 = qJD(1) * t213;
t209 = qJD(4) * t123;
t193 = t124 * t209;
t203 = qJD(2) * qJD(4);
t211 = qJD(2) * t128;
t79 = t127 * t203 + (t127 * t211 - t193) * qJD(1);
t40 = qJD(5) * t73 + t122 * t79 - t126 * t188;
t306 = -t40 / 0.2e1;
t379 = Ifges(6,4) * t306;
t300 = -t69 / 0.2e1;
t302 = -t46 / 0.2e1;
t304 = -t45 / 0.2e1;
t332 = t47 * mrSges(7,1) - t48 * mrSges(7,2);
t378 = Ifges(7,5) * t302 + Ifges(7,6) * t304 + Ifges(7,3) * t300 - t332;
t348 = t97 * mrSges(6,1);
t377 = t332 + t348 - t249 / 0.2e1 + t382;
t184 = qJD(6) * t126 - qJD(4);
t206 = qJD(5) * t122;
t222 = t123 * t125;
t208 = qJD(4) * t126;
t183 = -qJD(6) + t208;
t342 = t183 * t127;
t376 = t121 * (-t123 * t206 + t342) + t184 * t222;
t106 = (mrSges(4,1) * t128 + mrSges(4,3) * t124) * qJD(1);
t166 = t126 * qJD(3) - t101 * t122;
t236 = t122 * t166;
t202 = t123 * t236;
t221 = t123 * t126;
t232 = t127 * t97;
t139 = m(5) * (-t101 * t123 + t232) + m(6) * (-t221 * t82 + t202 + t232) + t106 + m(7) * t202;
t285 = mrSges(6,3) * t73;
t261 = -mrSges(6,1) * t103 - mrSges(7,1) * t45 + mrSges(7,2) * t46 + t285;
t53 = -mrSges(6,2) * t103 + t72 * mrSges(6,3);
t83 = -mrSges(5,2) * t116 + t104 * mrSges(5,3);
t154 = t122 * t261 + t126 * t53 + t83;
t256 = mrSges(5,3) * t105;
t260 = -mrSges(5,1) * t116 - mrSges(6,1) * t72 + mrSges(6,2) * t73 + t256;
t375 = -t123 * t154 + t127 * t260 + t139;
t207 = qJD(4) * t127;
t192 = t124 * t207;
t153 = t123 * t211 + t192;
t80 = qJD(1) * t153 + t123 * t203;
t294 = t80 / 0.2e1;
t39 = qJD(5) * t72 + t122 * t188 + t126 * t79;
t307 = t39 / 0.2e1;
t374 = Ifges(6,1) * t307 + Ifges(6,5) * t294;
t373 = t116 * Ifges(5,6) / 0.2e1;
t220 = t125 * t126;
t100 = -t121 * t127 + t123 * t220;
t151 = t128 * t100;
t219 = t126 * t127;
t159 = t121 * t123 + t125 * t219;
t165 = t184 * t127;
t217 = qJ(3) * qJD(1);
t204 = qJD(5) * t127;
t337 = t122 * t204 + t123 * t183;
t372 = t159 * qJD(3) + (-t121 * t165 - t125 * t337) * qJ(3) + t151 * t217;
t157 = t121 * t221 + t125 * t127;
t150 = t128 * t157;
t158 = -t121 * t219 + t222;
t371 = t158 * qJD(3) + (t121 * t337 - t125 * t165) * qJ(3) - t150 * t217;
t365 = qJD(5) * t166;
t210 = qJD(3) * t124;
t156 = -qJ(3) * t211 - t210;
t55 = (-qJ(3) * t209 + qJD(3) * t127) * qJD(2) + (-qJ(3) * t192 + t123 * t156) * qJD(1);
t32 = t126 * t55 + t365;
t56 = -qJD(4) * t182 - t156 * t215 + (qJ(3) * t207 + qJD(3) * t123) * qJD(2);
t7 = qJD(6) * t47 + t121 * t56 + t125 * t32;
t8 = -qJD(6) * t48 - t121 * t32 + t125 * t56;
t181 = t8 * mrSges(7,1) - t7 * mrSges(7,2);
t370 = t56 * mrSges(6,1) + t181;
t118 = Ifges(3,4) * t214;
t248 = t103 * Ifges(6,3);
t269 = t73 * Ifges(6,5);
t270 = t72 * Ifges(6,6);
t29 = t248 + t269 + t270;
t255 = Ifges(5,4) * t105;
t361 = t373 + t104 * t381 + t255 / 0.2e1;
t330 = -t166 * mrSges(6,1) + t82 * mrSges(6,2) + t361;
t138 = t29 / 0.2e1 + t270 / 0.2e1 + t269 / 0.2e1 + t248 / 0.2e1 - t101 * mrSges(5,3) - t330 - t361;
t344 = qJD(3) * mrSges(5,1);
t133 = (-t138 - t344) * t123;
t102 = Ifges(5,4) * t104;
t245 = t105 * Ifges(5,1);
t242 = t116 * Ifges(5,5);
t356 = -t242 / 0.2e1;
t61 = t102 + t242 + t245;
t145 = -t61 / 0.2e1 - t97 * mrSges(5,3) - t102 / 0.2e1 - t245 / 0.2e1 + t356;
t343 = qJD(3) * mrSges(5,2);
t142 = -t145 + t343;
t346 = qJD(2) / 0.2e1;
t368 = t142 * t127 + (t124 * Ifges(4,1) - Ifges(4,5) * t128) * qJD(1) / 0.2e1 + Ifges(3,1) * t216 / 0.2e1 + t118 / 0.2e1 + (Ifges(4,4) + Ifges(3,5)) * t346 - t133;
t250 = Ifges(6,5) * t103;
t284 = Ifges(6,1) * t73;
t68 = Ifges(6,4) * t72;
t31 = t250 + t68 + t284;
t347 = t97 * mrSges(6,2);
t367 = t31 / 0.2e1 + t347 + t250 / 0.2e1;
t15 = -qJD(6) * t46 - t121 * t39 + t125 * t80;
t12 = Ifges(7,6) * t15;
t14 = qJD(6) * t45 + t121 * t80 + t125 * t39;
t13 = Ifges(7,5) * t14;
t1 = Ifges(7,3) * t40 + t12 + t13;
t264 = t80 * Ifges(6,6);
t275 = t39 * Ifges(6,4);
t10 = -t40 * Ifges(6,2) + t264 + t275;
t364 = -t10 / 0.2e1 + t1 / 0.2e1 - t32 * mrSges(6,3);
t205 = qJD(5) * t126;
t33 = qJD(5) * t82 + t122 * t55;
t238 = t122 * t33;
t160 = -t166 * t205 + t238;
t229 = qJD(4) * t97;
t233 = t126 * t32;
t26 = -mrSges(6,2) * t80 - mrSges(6,3) * t40;
t288 = mrSges(6,1) * t80 + mrSges(7,1) * t15 - mrSges(7,2) * t14 - mrSges(6,3) * t39;
t362 = qJD(5) * (-t122 * t53 + t126 * t261) + m(5) * (t55 + t229) - m(6) * (t206 * t82 - t160 - t229 - t233) + m(7) * t160 + t260 * qJD(4) - t288 * t122 - mrSges(5,2) * t188 - mrSges(5,3) * t80 + t126 * t26;
t360 = t123 * (t121 * t184 + t125 * t206) - t125 * t342;
t359 = Ifges(5,5) / 0.2e1;
t358 = -Ifges(3,6) / 0.2e1;
t357 = -Ifges(5,6) / 0.2e1;
t316 = t14 / 0.2e1;
t315 = t15 / 0.2e1;
t305 = t40 / 0.2e1;
t352 = t33 * mrSges(6,3);
t350 = t56 * mrSges(6,2);
t349 = t166 * mrSges(6,3);
t176 = mrSges(7,1) * t121 + mrSges(7,2) * t125;
t345 = (mrSges(6,3) + t176) * t166;
t277 = t32 * mrSges(6,2);
t37 = Ifges(6,6) * t40;
t38 = Ifges(6,5) * t39;
t9 = Ifges(6,3) * t80 - t37 + t38;
t339 = (Ifges(6,3) / 0.2e1 + t381) * t80 + Ifges(5,2) * t294 + t188 * t357 - t79 * Ifges(5,4) + t9 / 0.2e1 + t38 / 0.2e1 - t37 / 0.2e1 - t277 - t33 * mrSges(6,1) - t55 * mrSges(5,3);
t331 = -t121 * t8 + t125 * t7;
t317 = t374 + t379;
t293 = -t103 / 0.2e1;
t329 = Ifges(6,5) * t293 - t347;
t327 = t56 * mrSges(5,1) + t55 * mrSges(5,2) - Ifges(5,5) * t79 + Ifges(5,6) * t80;
t325 = t317 + t350 + t374;
t323 = -Ifges(6,6) * t293 - t348 + t378;
t322 = -Ifges(7,5) * t316 + Ifges(6,2) * t306 + Ifges(6,6) * t294 - Ifges(7,6) * t315 - Ifges(7,3) * t305 - t370;
t2 = Ifges(7,4) * t14 + Ifges(7,2) * t15 + Ifges(7,6) * t40;
t320 = t2 / 0.2e1;
t319 = Ifges(7,1) * t316 + Ifges(7,4) * t315 + Ifges(7,5) * t305;
t282 = Ifges(7,4) * t46;
t18 = Ifges(7,2) * t45 + Ifges(7,6) * t69 + t282;
t312 = -t18 / 0.2e1;
t44 = Ifges(7,4) * t45;
t19 = Ifges(7,1) * t46 + Ifges(7,5) * t69 + t44;
t311 = -t19 / 0.2e1;
t303 = t45 / 0.2e1;
t301 = t46 / 0.2e1;
t299 = t69 / 0.2e1;
t298 = -t72 / 0.2e1;
t297 = t72 / 0.2e1;
t296 = -t73 / 0.2e1;
t295 = t73 / 0.2e1;
t289 = t127 / 0.2e1;
t287 = m(4) * qJ(3);
t286 = m(4) * qJ(3) ^ 2;
t194 = t123 * t214;
t218 = t126 * t128;
t88 = (t122 * t124 + t127 * t218) * qJD(1);
t67 = t121 * t194 + t125 * t88;
t259 = -t360 - t67;
t66 = -t121 * t88 + t125 * t194;
t258 = -t376 - t66;
t257 = m(4) * qJD(3);
t254 = Ifges(6,4) * t122;
t253 = Ifges(6,4) * t126;
t252 = Ifges(7,4) * t121;
t251 = Ifges(7,4) * t125;
t243 = t105 * t97;
t239 = t121 * t18;
t235 = t123 * t97;
t234 = t125 * t19;
t231 = qJD(2) * mrSges(4,3);
t230 = qJD(3) * mrSges(4,2);
t228 = t104 * t122;
t227 = t121 * t122;
t226 = t121 * t126;
t225 = t122 * t125;
t224 = t122 * t127;
t223 = t123 * t124;
t201 = t166 * t224;
t200 = -0.3e1 / 0.2e1 * Ifges(4,5) + 0.3e1 / 0.2e1 * Ifges(3,4);
t199 = Ifges(4,6) / 0.2e1 + t358;
t196 = qJD(4) * t236;
t180 = (Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1) * qJD(2);
t52 = (-qJD(5) + t212) * t218 + (-t123 * t208 + (qJD(2) - t204) * t122) * t124;
t179 = qJD(6) * t223 + t52;
t177 = mrSges(7,1) * t125 - mrSges(7,2) * t121;
t175 = Ifges(7,1) * t125 - t252;
t174 = Ifges(7,1) * t121 + t251;
t173 = -Ifges(7,2) * t121 + t251;
t172 = Ifges(7,2) * t125 + t252;
t171 = Ifges(7,5) * t125 - Ifges(7,6) * t121;
t170 = Ifges(7,5) * t121 + Ifges(7,6) * t125;
t168 = t121 * t47 - t125 * t48;
t167 = t126 * t82 - t236;
t27 = -mrSges(7,2) * t69 + mrSges(7,3) * t45;
t28 = mrSges(7,1) * t69 - mrSges(7,3) * t46;
t161 = -t121 * t28 + t125 * t27 + t53;
t152 = t124 * t100;
t149 = t157 * t124;
t148 = t56 * mrSges(5,3) + t79 * Ifges(5,1) - t80 * Ifges(5,4) + t188 * t359;
t99 = -t122 * t128 + t124 * t219;
t147 = -qJD(6) * t99 + t153;
t146 = t284 / 0.2e1 + t68 / 0.2e1 + t367;
t143 = t146 - t349;
t117 = Ifges(4,5) * t216;
t141 = t116 * Ifges(5,3) + t105 * Ifges(5,5) + t104 * Ifges(5,6) + Ifges(4,6) * t346 - Ifges(4,3) * t214 / 0.2e1 + t117 / 0.2e1 + qJD(2) * t358 - (Ifges(3,4) * t124 + t128 * Ifges(3,2)) * qJD(1) / 0.2e1 - t101 * mrSges(5,2) - t97 * mrSges(5,1);
t137 = -t280 / 0.2e1 - t283 / 0.2e1 + t279 / 0.2e1 + t281 / 0.2e1 + t278 / 0.2e1 + t377;
t135 = mrSges(6,1) * t40 + mrSges(6,2) * t39 - mrSges(5,1) * t188 + mrSges(5,3) * t79 - t154 * qJD(4) + m(7) * t196 + m(6) * (-t208 * t82 + t196 + t56) + m(5) * (-qJD(4) * t101 + t56);
t134 = t137 * t122;
t132 = -t380 - t239 / 0.2e1 + t234 / 0.2e1 + t173 * t303 + t175 * t301 + t171 * t299;
t130 = qJD(1) ^ 2;
t129 = qJD(2) ^ 2;
t120 = t124 ^ 2;
t115 = Ifges(5,3) * t188;
t109 = mrSges(4,2) * t214 + t231;
t90 = t159 * qJ(3);
t89 = t158 * qJ(3);
t86 = qJ(3) * t152;
t85 = qJ(3) * t149;
t71 = t121 * t223 + t125 * t99;
t70 = -t121 * t99 + t124 * t222;
t63 = t104 * t220 + t105 * t121;
t62 = -t104 * t226 + t105 * t125;
t58 = t101 * t121 - t220 * t97;
t57 = t101 * t125 + t226 * t97;
t24 = t121 * t147 + t125 * t179;
t23 = -t121 * t179 + t125 * t147;
t21 = qJD(3) * t149 + (qJD(2) * t150 + t124 * t376) * qJ(3);
t20 = -qJD(3) * t152 + (-qJD(2) * t151 + t124 * t360) * qJ(3);
t6 = -mrSges(7,2) * t40 + mrSges(7,3) * t15;
t5 = mrSges(7,1) * t40 - mrSges(7,3) * t14;
t3 = [t120 * t217 * t257 + m(7) * (t20 * t48 + t21 * t47 - t7 * t86 + t8 * t85) + (Ifges(7,1) * t71 + Ifges(7,4) * t70) * t316 + (Ifges(7,1) * t24 + Ifges(7,4) * t23) * t301 + (t23 * t48 - t24 * t47 + t7 * t70 - t71 * t8) * mrSges(7,3) + (Ifges(7,5) * t71 + Ifges(7,6) * t70) * t305 + (Ifges(7,5) * t24 + Ifges(7,6) * t23) * t299 + (Ifges(7,4) * t71 + Ifges(7,2) * t70) * t315 + (Ifges(7,4) * t24 + Ifges(7,2) * t23) * t303 - t166 * (-mrSges(7,1) * t23 + mrSges(7,2) * t24) + (t148 * t127 + t339 * t123 + (t123 * t145 + t127 * t138) * qJD(4) + (t199 * qJD(2) + ((Ifges(5,5) * t289 + t123 * t357 - t200) * t124 + (-Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + t286 + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t128) * qJD(1) + t141) * qJD(2) + (t106 + (mrSges(5,1) * qJD(4) + t260) * t127 + (-mrSges(5,2) * qJD(4) - t154) * t123 + t139) * qJD(3) + (-t129 * mrSges(4,2) + (m(4) * t210 + (-0.2e1 * mrSges(4,1) * t124 + (0.3e1 * mrSges(4,3) + t287) * t128) * qJD(2)) * qJD(1) + t135 * t127 - t362 * t123) * qJ(3)) * t124 + t71 * t319 + t70 * t320 + t23 * t18 / 0.2e1 + t24 * t19 / 0.2e1 + t20 * t27 + t21 * t28 + t33 * (-mrSges(7,1) * t70 + mrSges(7,2) * t71) + (-Ifges(6,4) * t307 - t322 + t364) * (t124 * t224 + t218) + t85 * t5 - t86 * t6 + (Ifges(6,1) * t295 + Ifges(6,4) * t297 - t349 + t367) * t52 + (-t115 / 0.2e1 + (0.2e1 * t230 + t180 + t200 * t214 + (mrSges(4,1) * t214 + t375) * qJ(3) + t368) * qJD(2) + t327) * t128 + (-Ifges(6,4) * t295 + Ifges(7,5) * t301 - Ifges(6,2) * t297 + Ifges(7,6) * t303 + Ifges(7,3) * t299 + t377) * (-t122 * t193 - t128 * t206 - t126 * t213 + (t122 * t211 + t124 * t205) * t127) + (t325 + t352 + t379) * t99; (t349 + Ifges(6,1) * t296 + Ifges(6,4) * t298 - t31 / 0.2e1 + t329) * t88 - t339 * t127 + (-t360 / 0.2e1 - t67 / 0.2e1) * t19 + (-Ifges(7,5) * t360 - Ifges(7,6) * t376) * t299 + (-Ifges(7,1) * t360 - Ifges(7,4) * t376) * t301 + (-Ifges(7,4) * t360 - Ifges(7,2) * t376) * t303 + (-t376 / 0.2e1 - t66 / 0.2e1) * t18 + (t231 + t109 + t260 * t123 + t154 * t127 + m(5) * (t101 * t127 + t235) + m(6) * (t219 * t82 - t201 + t235)) * qJD(3) + (Ifges(7,4) * t67 + Ifges(7,2) * t66) * t304 + (Ifges(7,4) * t100 - Ifges(7,2) * t157) * t315 + (Ifges(7,1) * t100 - Ifges(7,4) * t157) * t316 + t33 * (mrSges(7,1) * t157 + mrSges(7,2) * t100) - t157 * t320 + (Ifges(7,1) * t67 + Ifges(7,4) * t66) * t302 - (-mrSges(7,1) * t258 + mrSges(7,2) * t259) * t166 + (-t100 * t8 - t157 * t7 + t258 * t48 - t259 * t47) * mrSges(7,3) + (Ifges(7,5) * t100 - Ifges(7,6) * t157) * t305 + (-t133 + (t126 * t143 + t134 + t142) * t127) * qJD(4) + t100 * t319 + (-Ifges(6,4) * t296 - Ifges(6,2) * t298 + t323 - t382) * (-t126 * t216 + t214 * t224) - t124 * t130 * t128 * t286 + (0.2e1 * qJD(2) * t257 + t135 * t123 + t127 * t362) * qJ(3) + ((0.2e1 * t317 + t350 + t352) * t126 + (t13 / 0.2e1 + t12 / 0.2e1 - t275 / 0.2e1 - t264 / 0.2e1 + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t40 + t364 + t370) * t122 + (-t122 * t143 + t126 * t137) * qJD(5) + t148) * t123 + t89 * t5 + t90 * t6 + t371 * t28 + t372 * t27 + (-t201 * qJD(3) + t371 * t47 + t372 * t48 + t7 * t90 + t8 * t89) * m(7) + ((-t117 / 0.2e1 + (Ifges(3,4) / 0.2e1 + qJ(3) * mrSges(4,1)) * t216 + (Ifges(5,6) * t289 + t123 * t359 + t199) * qJD(2) - t141) * t124 + (-t230 - t118 / 0.2e1 + Ifges(4,5) * t214 / 0.2e1 + t180 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t216 + (-mrSges(4,3) * t216 - t375) * qJ(3) - t368) * t128) * qJD(1) + (Ifges(7,5) * t67 + Ifges(7,6) * t66) * t300; t80 * mrSges(5,1) + t79 * mrSges(5,2) - qJD(2) * t109 - t104 * t83 - t63 * t27 - t62 * t28 - t260 * t105 + (mrSges(4,2) * t211 - t106 * t124) * qJD(1) + (qJD(5) * t161 - t104 * t53 + t288) * t126 + (-t121 * t5 + t125 * t6 + t26 + (-t121 * t27 - t125 * t28) * qJD(6) + t103 * t261) * t122 - m(5) * (t101 * t104 + t243) - (t120 * t130 + t129) * t287 + ((-qJD(5) * t168 - t33) * t126 + (-qJD(6) * t328 + t331 - t365) * t122 + t228 * t166 - t47 * t62 - t48 * t63) * m(7) + (t103 * t167 + t122 * t32 - t126 * t33 - t243) * m(6); (t256 - t260) * t101 + ((t132 + t146 - t345) * qJD(5) + t10 / 0.2e1 + t322) * t126 + t134 * qJD(5) - t2 * t227 / 0.2e1 - m(7) * (t47 * t57 + t48 * t58) + t176 * t238 + t253 * t307 + (t238 + t233) * mrSges(6,3) + t115 + (-t225 * t8 - t227 * t7) * mrSges(7,3) - (t104 * t31 + t1) * t126 / 0.2e1 - (Ifges(5,1) * t104 - t255 + t29) * t105 / 0.2e1 - (-Ifges(5,2) * t105 + t102 + t61) * t104 / 0.2e1 + t254 * t306 - t327 + ((mrSges(7,3) * t168 + t121 * t311 + t125 * t312 - t166 * t177 + t170 * t300 + t172 * t304 + t174 * t302) * qJD(6) + t171 * t305 + t173 * t315 + t175 * t316 + t325) * t122 + ((Ifges(6,5) * t126 - Ifges(6,6) * t122) * t293 + (Ifges(6,1) * t126 - t254) * t296 + (-Ifges(6,2) * t122 + t253) * t298 + (-mrSges(6,1) * t122 - mrSges(6,2) * t126 - mrSges(5,3)) * t97 + (t122 * t82 + t126 * t166) * mrSges(6,3) - t343 + t356) * t104 + (mrSges(7,2) * t166 + mrSges(7,3) * t47 + Ifges(7,1) * t302 + Ifges(7,4) * t304 + Ifges(7,5) * t300 + t311) * t63 + (-mrSges(7,1) * t166 - mrSges(7,3) * t48 + Ifges(7,4) * t302 + Ifges(7,2) * t304 + Ifges(7,6) * t300 + t312) * t62 + (-m(6) * (t101 - t167) - m(7) * t236 + t154) * t97 + t225 * t319 - t57 * t28 - t58 * t27 + (Ifges(6,5) * t296 + Ifges(6,6) * t298 + Ifges(6,3) * t293 + t330 - t344 + t373) * t105 + (t369 + t378) * t228; (-t261 + t285) * t82 + t170 * t305 + t172 * t315 + t174 * t316 + t331 * mrSges(7,3) - t277 + t9 + t30 * t295 + t239 * t297 + (-t166 * t176 + t132) * qJD(6) + (-mrSges(6,1) - t177) * t33 - (-m(7) * (t168 + t82) + t161) * t166 + t121 * t319 + t125 * t320 + (t171 * t300 + t173 * t304 + t175 * t302 + t329 + t345 + t380) * t72 + t323 * t73 + (Ifges(6,1) * t72 + t17 - t283) * t296 + (-Ifges(6,2) * t73 + t234 + t31 + t68) * t298; t166 * (mrSges(7,1) * t46 + mrSges(7,2) * t45) + (Ifges(7,1) * t45 - t282) * t302 + t18 * t301 + (Ifges(7,5) * t45 - Ifges(7,6) * t46) * t300 - t47 * t27 + t48 * t28 + (t45 * t47 + t46 * t48) * mrSges(7,3) + t181 + t1 + (-Ifges(7,2) * t46 + t19 + t44) * t304;];
tauc  = t3(:);

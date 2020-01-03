% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:48
% EndTime: 2019-12-31 22:04:17
% DurationCPUTime: 12.66s
% Computational Cost: add. (5516->577), mult. (13837->801), div. (0->0), fcn. (8928->6), ass. (0->241)
t321 = Ifges(5,1) + Ifges(6,1);
t320 = Ifges(6,4) + Ifges(5,5);
t319 = Ifges(6,5) - Ifges(5,4);
t318 = Ifges(6,2) + Ifges(5,3);
t317 = -Ifges(5,6) + Ifges(6,6);
t193 = sin(qJ(4));
t194 = sin(qJ(3));
t196 = cos(qJ(4));
t197 = cos(qJ(3));
t153 = t193 * t194 - t196 * t197;
t311 = qJD(3) + qJD(4);
t109 = t311 * t153;
t198 = cos(qJ(2));
t204 = t153 * t198;
t129 = qJD(1) * t204;
t328 = t109 - t129;
t154 = t193 * t197 + t194 * t196;
t110 = t311 * t154;
t256 = qJD(1) * t198;
t128 = t154 * t256;
t260 = t110 - t128;
t181 = qJD(3) - t256;
t170 = qJD(4) + t181;
t195 = sin(qJ(2));
t257 = qJD(1) * t195;
t238 = t197 * t257;
t152 = qJD(2) * t194 + t238;
t239 = t194 * t257;
t254 = qJD(2) * t197;
t205 = t239 - t254;
t201 = t196 * t152 - t193 * t205;
t101 = t193 * t152 + t196 * t205;
t271 = Ifges(6,5) * t101;
t95 = Ifges(5,4) * t101;
t314 = t320 * t170 + t321 * t201 + t271 - t95;
t327 = t314 / 0.2e1;
t251 = qJD(3) * t195;
t253 = qJD(2) * t198;
t203 = -t194 * t251 + t197 * t253;
t246 = qJD(2) * qJD(3);
t118 = qJD(1) * t203 + t197 * t246;
t237 = t194 * t253;
t250 = qJD(3) * t197;
t202 = t195 * t250 + t237;
t119 = -qJD(1) * t202 - t194 * t246;
t39 = qJD(4) * t201 + t193 * t118 - t196 * t119;
t304 = t39 / 0.2e1;
t326 = Ifges(6,3) * t304;
t56 = pkin(4) * t201 + qJ(5) * t101;
t325 = -Ifges(5,6) / 0.2e1;
t38 = -qJD(4) * t101 + t196 * t118 + t193 * t119;
t306 = t38 / 0.2e1;
t295 = t118 / 0.2e1;
t294 = t119 / 0.2e1;
t247 = qJD(1) * qJD(2);
t231 = t195 * t247;
t323 = t231 / 0.2e1;
t316 = t320 * t231 + t319 * t39 + t321 * t38;
t190 = pkin(6) * t256;
t283 = pkin(3) * t194;
t148 = t256 * t283 + t190;
t252 = qJD(3) * t194;
t315 = pkin(3) * t252 + t260 * pkin(4) + t328 * qJ(5) - qJD(5) * t154 - t148;
t134 = t154 * t195;
t160 = -pkin(2) * t198 - t195 * pkin(7) - pkin(1);
t145 = t160 * qJD(1);
t167 = qJD(2) * pkin(7) + t190;
t107 = t194 * t145 + t197 * t167;
t81 = -pkin(8) * t205 + t107;
t268 = t193 * t81;
t106 = t197 * t145 - t167 * t194;
t80 = -pkin(8) * t152 + t106;
t74 = pkin(3) * t181 + t80;
t21 = t196 * t74 - t268;
t313 = qJD(5) - t21;
t262 = t197 * t198;
t183 = pkin(6) * t262;
t127 = t194 * t160 + t183;
t312 = t318 * t231 + t317 * t39 + t320 * t38;
t189 = pkin(6) * t257;
t166 = -qJD(2) * pkin(2) + t189;
t210 = t106 * t197 + t107 * t194;
t270 = Ifges(4,6) * t194;
t215 = Ifges(4,5) * t197 - t270;
t276 = Ifges(4,4) * t194;
t219 = Ifges(4,1) * t197 - t276;
t220 = mrSges(4,1) * t194 + mrSges(4,2) * t197;
t286 = t197 / 0.2e1;
t287 = t181 / 0.2e1;
t291 = t152 / 0.2e1;
t277 = Ifges(4,4) * t152;
t88 = -Ifges(4,2) * t205 + Ifges(4,6) * t181 + t277;
t149 = Ifges(4,4) * t205;
t89 = t152 * Ifges(4,1) + t181 * Ifges(4,5) - t149;
t310 = -t210 * mrSges(4,3) + t166 * t220 + t219 * t291 + t215 * t287 - t194 * t88 / 0.2e1 + t89 * t286;
t267 = t196 * t81;
t22 = t193 * t74 + t267;
t221 = pkin(2) * t195 - pkin(7) * t198;
t158 = t221 * qJD(2);
t146 = qJD(1) * t158;
t225 = pkin(6) * t231;
t65 = -qJD(3) * t107 + t197 * t146 + t194 * t225;
t30 = pkin(3) * t231 - pkin(8) * t118 + t65;
t64 = t145 * t250 + t194 * t146 - t167 * t252 - t197 * t225;
t40 = pkin(8) * t119 + t64;
t6 = -qJD(4) * t22 - t193 * t40 + t196 * t30;
t151 = t197 * t160;
t263 = t195 * t197;
t282 = pkin(6) * t194;
t105 = -pkin(8) * t263 + t151 + (-pkin(3) - t282) * t198;
t265 = t194 * t195;
t112 = -pkin(8) * t265 + t127;
t266 = t193 * t105 + t196 * t112;
t208 = pkin(3) * t195 - pkin(8) * t262;
t255 = qJD(2) * t195;
t258 = t197 * t158 + t255 * t282;
t59 = t208 * qJD(2) + (-t183 + (pkin(8) * t195 - t160) * t194) * qJD(3) + t258;
t78 = t194 * t158 + t160 * t250 + (-t195 * t254 - t198 * t252) * pkin(6);
t63 = -pkin(8) * t202 + t78;
t11 = -qJD(4) * t266 - t193 * t63 + t196 * t59;
t124 = pkin(3) * t205 + t166;
t20 = qJ(5) * t170 + t22;
t31 = t101 * pkin(4) - qJ(5) * t201 + t124;
t94 = Ifges(6,5) * t201;
t45 = Ifges(6,6) * t170 + Ifges(6,3) * t101 + t94;
t274 = Ifges(5,4) * t201;
t48 = -Ifges(5,2) * t101 + Ifges(5,6) * t170 + t274;
t309 = t20 * mrSges(6,2) + t22 * mrSges(5,3) - t45 / 0.2e1 + t48 / 0.2e1 - t124 * mrSges(5,1) - t31 * mrSges(6,1);
t308 = Ifges(6,5) * t306 + Ifges(6,6) * t323 + t326;
t307 = -t38 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t304 + t231 * t325;
t305 = -t39 / 0.2e1;
t302 = Ifges(4,1) * t295 + Ifges(4,4) * t294 + Ifges(4,5) * t323;
t301 = -pkin(8) - pkin(7);
t300 = -t101 / 0.2e1;
t299 = t101 / 0.2e1;
t297 = -t201 / 0.2e1;
t296 = t201 / 0.2e1;
t292 = -t152 / 0.2e1;
t289 = -t170 / 0.2e1;
t288 = t170 / 0.2e1;
t284 = m(5) * t124;
t82 = -mrSges(6,2) * t101 + mrSges(6,3) * t170;
t279 = mrSges(5,3) * t101;
t83 = -mrSges(5,2) * t170 - t279;
t281 = t82 + t83;
t278 = mrSges(5,3) * t201;
t84 = mrSges(5,1) * t170 - t278;
t85 = -mrSges(6,1) * t170 + mrSges(6,2) * t201;
t280 = t84 - t85;
t156 = t221 * qJD(1);
t140 = t194 * t156;
t264 = t194 * t198;
t104 = t140 + (-pkin(6) * t263 - pkin(8) * t264) * qJD(1);
t120 = pkin(6) * t239 + t197 * t156;
t90 = qJD(1) * t208 + t120;
t52 = t196 * t104 + t193 * t90;
t275 = Ifges(4,4) * t197;
t273 = Ifges(4,5) * t152;
t272 = Ifges(4,5) * t194;
t269 = Ifges(4,3) * t181;
t259 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t205 + t152 * mrSges(4,2) + mrSges(3,3) * t257;
t159 = pkin(3) * t265 + t195 * pkin(6);
t249 = qJD(4) * t193;
t248 = qJD(4) * t196;
t191 = pkin(6) * t253;
t243 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t242 = Ifges(6,6) / 0.2e1 + t325;
t241 = Ifges(4,5) * t118 + Ifges(4,6) * t119 + Ifges(4,3) * t231;
t125 = pkin(3) * t202 + t191;
t187 = -pkin(3) * t197 - pkin(2);
t240 = qJD(3) * t301;
t230 = t198 * t247;
t228 = t253 / 0.2e1;
t227 = -t251 / 0.2e1;
t99 = -pkin(3) * t119 + pkin(6) * t230;
t224 = t197 * t240;
t218 = Ifges(4,1) * t194 + t275;
t217 = -Ifges(4,2) * t194 + t275;
t216 = Ifges(4,2) * t197 + t276;
t212 = -t194 * t65 + t197 * t64;
t51 = -t104 * t193 + t196 * t90;
t60 = t196 * t105 - t193 * t112;
t168 = t301 * t194;
t169 = t301 * t197;
t209 = t196 * t168 + t169 * t193;
t117 = t168 * t193 - t169 * t196;
t25 = -mrSges(6,1) * t231 + t38 * mrSges(6,2);
t5 = t193 * t30 + t196 * t40 + t74 * t248 - t249 * t81;
t10 = t105 * t248 - t112 * t249 + t193 * t59 + t196 * t63;
t207 = t194 * t217;
t206 = t197 * t217;
t2 = qJ(5) * t231 + qJD(5) * t170 + t5;
t3 = -pkin(4) * t231 - t6;
t200 = t6 * mrSges(5,1) - t3 * mrSges(6,1) - t5 * mrSges(5,2) + t2 * mrSges(6,3) + t312;
t188 = Ifges(3,4) * t256;
t186 = -pkin(3) * t196 - pkin(4);
t184 = pkin(3) * t193 + qJ(5);
t176 = pkin(3) * t248 + qJD(5);
t164 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t256;
t157 = t194 * t240;
t139 = Ifges(3,1) * t257 + Ifges(3,5) * qJD(2) + t188;
t138 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t195 + t198 * Ifges(3,2)) * qJD(1);
t135 = t153 * t195;
t126 = -pkin(6) * t264 + t151;
t123 = mrSges(4,1) * t181 - mrSges(4,3) * t152;
t122 = -t181 * mrSges(4,2) - mrSges(4,3) * t205;
t121 = -pkin(6) * t238 + t140;
t98 = -mrSges(4,2) * t231 + mrSges(4,3) * t119;
t97 = mrSges(4,1) * t231 - mrSges(4,3) * t118;
t96 = pkin(4) * t153 - qJ(5) * t154 + t187;
t87 = -Ifges(4,6) * t205 + t269 + t273;
t79 = -qJD(3) * t127 + t258;
t76 = pkin(4) * t134 + qJ(5) * t135 + t159;
t75 = -mrSges(4,1) * t119 + mrSges(4,2) * t118;
t73 = qJD(4) * t117 + t157 * t193 - t196 * t224;
t72 = qJD(4) * t209 + t196 * t157 + t193 * t224;
t68 = t118 * Ifges(4,4) + t119 * Ifges(4,2) + Ifges(4,6) * t231;
t67 = -t249 * t265 + (t263 * t311 + t237) * t196 + t203 * t193;
t66 = -qJD(2) * t204 - t134 * t311;
t58 = mrSges(5,1) * t101 + mrSges(5,2) * t201;
t57 = mrSges(6,1) * t101 - mrSges(6,3) * t201;
t55 = pkin(4) * t198 - t60;
t54 = -qJ(5) * t198 + t266;
t47 = Ifges(6,4) * t201 + t170 * Ifges(6,2) + t101 * Ifges(6,6);
t46 = Ifges(5,5) * t201 - t101 * Ifges(5,6) + t170 * Ifges(5,3);
t44 = -pkin(4) * t257 - t51;
t43 = qJ(5) * t257 + t52;
t42 = pkin(3) * t152 + t56;
t28 = t196 * t80 - t268;
t27 = t193 * t80 + t267;
t26 = -mrSges(5,2) * t231 - mrSges(5,3) * t39;
t24 = mrSges(5,1) * t231 - mrSges(5,3) * t38;
t23 = -mrSges(6,2) * t39 + mrSges(6,3) * t231;
t19 = -pkin(4) * t170 + t313;
t18 = pkin(4) * t67 - qJ(5) * t66 + qJD(5) * t135 + t125;
t17 = mrSges(5,1) * t39 + mrSges(5,2) * t38;
t16 = mrSges(6,1) * t39 - mrSges(6,3) * t38;
t9 = -pkin(4) * t255 - t11;
t8 = qJ(5) * t255 - qJD(5) * t198 + t10;
t7 = pkin(4) * t39 - qJ(5) * t38 - qJD(5) * t201 + t99;
t1 = [((0.3e1 / 0.2e1 * t198 ^ 2 - 0.3e1 / 0.2e1 * t195 ^ 2) * Ifges(3,4) - 0.2e1 * (mrSges(3,1) * t195 + mrSges(3,2) * t198) * pkin(1)) * t247 + m(4) * (t106 * t79 + t107 * t78 + t65 * t126 + t64 * t127 + (t166 + t189) * t191) + (Ifges(5,4) * t66 + Ifges(5,6) * t255) * t300 - t205 * (-t216 * t251 + (Ifges(4,6) * t195 + t198 * t217) * qJD(2)) / 0.2e1 + (-t210 * t253 + (-t64 * t194 - t65 * t197 + (t106 * t194 - t107 * t197) * qJD(3)) * t195) * mrSges(4,3) - t138 * t255 / 0.2e1 + (-t107 * t255 + t166 * t203 + t64 * t198) * mrSges(4,2) + t19 * (-mrSges(6,1) * t255 + mrSges(6,2) * t66) + (-Ifges(5,4) * t135 - Ifges(5,6) * t198) * t305 + (-Ifges(6,5) * t135 - Ifges(6,6) * t198) * t304 + (t124 * t66 - t135 * t99 + t198 * t5 - t22 * t255) * mrSges(5,2) + t3 * (mrSges(6,1) * t198 - t135 * mrSges(6,2)) + t6 * (-mrSges(5,1) * t198 + t135 * mrSges(5,3)) + (t135 * t7 - t198 * t2 + t20 * t255 - t31 * t66) * mrSges(6,3) + t266 * t26 + m(5) * (t10 * t22 + t11 * t21 + t124 * t125 + t159 * t99 + t266 * t5 + t6 * t60) + (Ifges(6,5) * t66 + Ifges(6,6) * t255) * t299 + (t194 * t227 + t197 * t228) * t89 + (-Ifges(4,6) * t198 + t195 * t217) * t294 + (-Ifges(4,5) * t198 + t195 * t219) * t295 + t263 * t302 + ((-Ifges(4,6) * t197 - t272) * t251 + (Ifges(4,3) * t195 + t198 * t215) * qJD(2)) * t287 + (-t218 * t251 + (Ifges(4,5) * t195 + t198 * t219) * qJD(2)) * t291 + (t197 * t227 - t237 / 0.2e1) * t88 - t68 * t265 / 0.2e1 + (mrSges(5,1) * t99 + mrSges(6,1) * t7 - mrSges(6,2) * t2 - mrSges(5,3) * t5 - Ifges(5,2) * t305 + t306 * t319 + t307 + t308 + t326) * t134 + (t195 * t75 + (t259 * t198 + (t220 * t256 - t164) * t195) * qJD(2)) * pkin(6) + (t106 * t255 + t166 * t202 - t65 * t198) * mrSges(4,1) + t21 * (mrSges(5,1) * t255 - mrSges(5,3) * t66) + ((t195 * t215 - t320 * t135 + t317 * t134 + (-Ifges(4,3) - t318) * t198) * qJD(1) + t87 + t47 + t46) * t255 / 0.2e1 - (t241 + t312) * t198 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t195 * t230 + qJD(2) ^ 2 * (Ifges(3,5) * t198 - Ifges(3,6) * t195) / 0.2e1 - t316 * t135 / 0.2e1 + (-Ifges(5,2) * t300 + Ifges(6,3) * t299 + t317 * t288 + t319 * t296 - t309) * t67 + (t255 * t318 + t320 * t66) * t288 + (t255 * t320 + t321 * t66) * t296 + (-t135 * t321 - t198 * t320) * t306 + t139 * t228 + t54 * t23 + t55 * t25 + t18 * t57 + t60 * t24 + m(6) * (t18 * t31 + t19 * t9 + t2 * t54 + t20 * t8 + t3 * t55 + t7 * t76) + t76 * t16 + t8 * t82 + t10 * t83 + t11 * t84 + t9 * t85 + t78 * t122 + t79 * t123 + t125 * t58 + t126 * t97 + t127 * t98 + t159 * t17 + t66 * t327; (m(4) * t212 + (-m(4) * t210 - t194 * t122 - t197 * t123) * qJD(3) - t194 * t97 + t197 * t98) * pkin(7) - m(4) * (t106 * t120 + t107 * t121) + (t45 - t48) * (t110 / 0.2e1 - t128 / 0.2e1) + t281 * t72 + (mrSges(5,1) * t260 - mrSges(5,2) * t328) * t124 + (-t153 * t2 + t154 * t3 - t19 * t328 - t20 * t260) * mrSges(6,2) + (mrSges(6,1) * t260 + mrSges(6,3) * t328) * t31 + (-t153 * t5 - t154 * t6 + t21 * t328 - t22 * t260) * mrSges(5,3) - t280 * t73 + t216 * t294 + t218 * t295 + t194 * t302 + t68 * t286 + (t209 * t6 + t117 * t5 - t124 * t148 + t187 * t99 + (-t52 + t72) * t22 + (-t51 - t73) * t21) * m(5) + (-t209 * t3 + t117 * t2 + t7 * t96 + t315 * t31 + (-t43 + t72) * t20 + (-t44 + t73) * t19) * m(6) - (t25 - t24) * t209 + (-t109 * t320 + t110 * t317) * t288 + (-t109 * t321 + t110 * t319) * t296 + (-Ifges(5,4) * t109 - Ifges(6,5) * t129 - Ifges(5,2) * t110 + Ifges(6,3) * t128) * t300 + (-Ifges(5,4) * t129 - Ifges(6,5) * t109 - Ifges(5,2) * t128 + Ifges(6,3) * t110) * t299 + t314 * (-t109 / 0.2e1 + t129 / 0.2e1) + (t128 * t317 - t129 * t320) * t289 + (t128 * t319 - t129 * t321) * t297 + (Ifges(6,5) * t154 + Ifges(6,3) * t153) * t304 + (Ifges(5,4) * t154 - Ifges(5,2) * t153) * t305 + t153 * t307 + t153 * t308 + ((-t139 / 0.2e1 - t188 / 0.2e1 + pkin(1) * mrSges(3,2) * qJD(1) + (Ifges(3,5) / 0.2e1 - t206 / 0.2e1) * qJD(2) + (-m(4) * t166 + (-m(4) * pkin(2) - mrSges(4,1) * t197 + mrSges(4,2) * t194 - mrSges(3,1)) * qJD(2) - t259) * pkin(6) - t310) * t198 + (t138 / 0.2e1 - t87 / 0.2e1 - t47 / 0.2e1 - t46 / 0.2e1 - t106 * mrSges(4,1) + t107 * mrSges(4,2) + t19 * mrSges(6,1) - t21 * mrSges(5,1) - t20 * mrSges(6,3) + t22 * mrSges(5,2) + pkin(6) * t164 - t273 / 0.2e1 - t269 / 0.2e1 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t170 - t243 * t201 - t242 * t101 - qJD(3) * t207 / 0.2e1 + ((t270 / 0.2e1 + Ifges(3,4) / 0.2e1) * t195 + pkin(1) * mrSges(3,1)) * qJD(1) + (t207 / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t256 + (pkin(6) * mrSges(3,2) + t272 / 0.2e1 - Ifges(3,6) / 0.2e1 + t243 * t154 + t242 * t153) * qJD(2)) * t195) * qJD(1) + (qJD(2) * t206 / 0.2e1 + (t58 + t284) * t283 + t310) * qJD(3) + t315 * t57 + t316 * t154 / 0.2e1 + (t153 * t319 + t154 * t321) * t306 + t212 * mrSges(4,3) + (t23 + t26) * t117 - pkin(2) * t75 - t43 * t82 - t52 * t83 - t51 * t84 - t44 * t85 + t96 * t16 - t121 * t122 - t120 * t123 - t148 * t58 + t7 * (mrSges(6,1) * t153 - mrSges(6,3) * t154) + t99 * (mrSges(5,1) * t153 + mrSges(5,2) * t154) + t187 * t17; t200 - t281 * t28 + t280 * t27 + (-Ifges(4,1) * t205 - t277) * t292 + t88 * t291 + t241 - m(5) * (-t21 * t27 + t22 * t28) + (-t152 * t58 + t193 * t26 + t196 * t24 + (-t193 * t280 + t196 * t83) * qJD(4) + 0.2e1 * t284 * t292 + m(6) * t19 * t249 + m(5) * (t193 * t5 + t196 * t6 - t21 * t249 + t22 * t248)) * pkin(3) + (-t106 * t205 + t107 * t152) * mrSges(4,3) - t166 * (t152 * mrSges(4,1) - mrSges(4,2) * t205) - t181 * (-Ifges(4,5) * t205 - Ifges(4,6) * t152) / 0.2e1 + (t184 * t2 + t186 * t3 - t19 * t27 - t31 * t42 + (t176 - t28) * t20) * m(6) + (-Ifges(5,2) * t299 + Ifges(6,3) * t300 + t317 * t289 + t319 * t297 + t309) * t201 - t42 * t57 - t64 * mrSges(4,2) + t65 * mrSges(4,1) - t106 * t122 + t107 * t123 + t176 * t82 + t184 * t23 + t186 * t25 + (-Ifges(4,2) * t152 - t149 + t89) * t205 / 0.2e1 + (mrSges(5,2) * t124 + mrSges(6,2) * t19 - mrSges(5,3) * t21 - mrSges(6,3) * t31 - Ifges(5,4) * t299 - Ifges(6,5) * t300 - t320 * t289 - t321 * t297 + t327) * t101; (t278 + t280) * t22 + (-t279 - t281) * t21 + (t101 * t19 + t20 * t201) * mrSges(6,2) + t200 + qJ(5) * t23 - pkin(4) * t25 - t56 * t57 + qJD(5) * t82 - t31 * (mrSges(6,1) * t201 + mrSges(6,3) * t101) + (Ifges(6,3) * t201 - t271) * t300 + t48 * t296 - t124 * (mrSges(5,1) * t201 - mrSges(5,2) * t101) + (-t101 * t320 + t201 * t317) * t289 + (-pkin(4) * t3 + qJ(5) * t2 - t19 * t22 + t20 * t313 - t31 * t56) * m(6) + (-Ifges(5,2) * t201 + t314 - t95) * t299 + (-t101 * t321 - t274 + t45 + t94) * t297; t201 * t57 - t170 * t82 + 0.2e1 * (t3 / 0.2e1 + t31 * t296 + t20 * t289) * m(6) + t25;];
tauc = t1(:);

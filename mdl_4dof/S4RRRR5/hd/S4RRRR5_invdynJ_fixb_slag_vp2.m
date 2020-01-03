% Calculate vector of inverse dynamics joint torques for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:23
% EndTime: 2019-12-31 17:27:44
% DurationCPUTime: 11.89s
% Computational Cost: add. (3439->488), mult. (7853->695), div. (0->0), fcn. (5044->10), ass. (0->238)
t169 = cos(qJ(2));
t234 = qJD(1) * t169;
t148 = qJD(3) - t234;
t144 = qJD(4) + t148;
t279 = t144 / 0.2e1;
t164 = sin(qJ(3));
t165 = sin(qJ(2));
t235 = qJD(1) * t165;
t217 = t164 * t235;
t168 = cos(qJ(3));
t232 = qJD(2) * t168;
t125 = -t217 + t232;
t215 = t168 * t235;
t126 = qJD(2) * t164 + t215;
t163 = sin(qJ(4));
t167 = cos(qJ(4));
t64 = t125 * t163 + t126 * t167;
t285 = t64 / 0.2e1;
t205 = t167 * t125 - t126 * t163;
t287 = t205 / 0.2e1;
t334 = Ifges(5,5) * t285 + Ifges(5,6) * t287 + Ifges(5,3) * t279;
t333 = t125 / 0.2e1;
t281 = t126 / 0.2e1;
t332 = t148 / 0.2e1;
t331 = -qJD(1) / 0.2e1;
t330 = mrSges(4,3) + mrSges(5,3);
t329 = Ifges(3,2) / 0.2e1;
t227 = qJD(1) * qJD(2);
t134 = qJDD(1) * t165 + t169 * t227;
t328 = t134 / 0.2e1;
t231 = qJD(2) * t169;
t220 = pkin(5) * t231;
t203 = pkin(2) * t165 - pkin(6) * t169;
t129 = t203 * qJD(1);
t109 = t164 * t129;
t171 = -pkin(7) - pkin(6);
t218 = qJD(3) * t171;
t241 = t165 * t168;
t243 = t164 * t169;
t327 = t164 * t218 - t109 - (-pkin(5) * t241 - pkin(7) * t243) * qJD(1);
t239 = t168 * t169;
t187 = pkin(3) * t165 - pkin(7) * t239;
t76 = pkin(5) * t217 + t168 * t129;
t326 = -qJD(1) * t187 + t168 * t218 - t76;
t151 = pkin(3) * t168 + pkin(2);
t162 = qJ(3) + qJ(4);
t157 = sin(t162);
t158 = cos(t162);
t200 = -mrSges(4,1) * t168 + mrSges(4,2) * t164;
t325 = m(4) * pkin(2) + m(5) * t151 + mrSges(5,1) * t158 - mrSges(5,2) * t157 - t200;
t324 = -m(4) * pkin(6) + m(5) * t171 - t330;
t201 = mrSges(3,1) * t165 + mrSges(3,2) * t169;
t262 = Ifges(3,4) * t165;
t323 = -pkin(1) * t201 + t165 * (Ifges(3,1) * t169 - t262) / 0.2e1;
t204 = pkin(2) * t169 + pkin(6) * t165;
t136 = -pkin(1) - t204;
t115 = t136 * qJD(1);
t156 = pkin(5) * t234;
t141 = qJD(2) * pkin(6) + t156;
t228 = qJD(3) * t168;
t230 = qJD(3) * t164;
t133 = qJDD(1) * t169 - t165 * t227;
t249 = qJDD(1) * pkin(1);
t69 = -pkin(2) * t133 - pkin(6) * t134 - t249;
t120 = t133 * pkin(5);
t99 = qJDD(2) * pkin(6) + t120;
t23 = t115 * t228 - t141 * t230 + t164 * t69 + t168 * t99;
t57 = -qJD(3) * t126 + qJDD(2) * t168 - t134 * t164;
t12 = pkin(7) * t57 + t23;
t68 = t115 * t164 + t141 * t168;
t47 = pkin(7) * t125 + t68;
t256 = t163 * t47;
t67 = t168 * t115 - t141 * t164;
t46 = -pkin(7) * t126 + t67;
t40 = pkin(3) * t148 + t46;
t17 = t167 * t40 - t256;
t122 = qJDD(3) - t133;
t24 = -t68 * qJD(3) - t164 * t99 + t168 * t69;
t56 = qJD(3) * t125 + qJDD(2) * t164 + t134 * t168;
t9 = pkin(3) * t122 - pkin(7) * t56 + t24;
t2 = qJD(4) * t17 + t12 * t167 + t163 * t9;
t252 = t167 * t47;
t18 = t163 * t40 + t252;
t3 = -qJD(4) * t18 - t12 * t163 + t167 * t9;
t322 = -t3 * mrSges(5,1) + t2 * mrSges(5,2);
t321 = -t24 * mrSges(4,1) + t23 * mrSges(4,2);
t114 = qJDD(4) + t122;
t284 = t114 / 0.2e1;
t16 = -qJD(4) * t64 - t163 * t56 + t167 * t57;
t297 = t16 / 0.2e1;
t15 = qJD(4) * t205 + t163 * t57 + t167 * t56;
t298 = t15 / 0.2e1;
t300 = Ifges(5,1) * t298 + Ifges(5,4) * t297 + Ifges(5,5) * t284;
t301 = Ifges(5,4) * t298 + Ifges(5,2) * t297 + Ifges(5,6) * t284;
t195 = t169 * Ifges(3,2) + t262;
t320 = t17 * mrSges(5,1) + t67 * mrSges(4,1) - Ifges(3,6) * qJD(2) / 0.2e1 + t195 * t331 + Ifges(4,5) * t281 + Ifges(4,6) * t333 + Ifges(4,3) * t332 - t18 * mrSges(5,2) - t68 * mrSges(4,2) + t334;
t299 = m(5) * pkin(3);
t290 = t56 / 0.2e1;
t289 = t57 / 0.2e1;
t283 = t122 / 0.2e1;
t318 = mrSges(2,2) - mrSges(3,3);
t142 = t171 * t164;
t143 = t171 * t168;
t74 = t142 * t167 + t143 * t163;
t317 = qJD(4) * t74 + t163 * t326 + t167 * t327;
t75 = t142 * t163 - t143 * t167;
t316 = -qJD(4) * t75 - t163 * t327 + t167 * t326;
t314 = mrSges(4,1) + t299;
t313 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t125 - mrSges(4,2) * t126 - mrSges(3,3) * t235;
t188 = t163 * t164 - t167 * t168;
t98 = t188 * t165;
t154 = Ifges(3,4) * t234;
t119 = Ifges(4,4) * t125;
t52 = t126 * Ifges(4,1) + t148 * Ifges(4,5) + t119;
t312 = Ifges(3,1) * t235 + Ifges(3,5) * qJD(2) + t168 * t52 + t154;
t216 = t164 * t234;
t311 = -t156 + (-t216 + t230) * pkin(3);
t149 = pkin(5) * t239;
t83 = t164 * t136 + t149;
t121 = t134 * pkin(5);
t310 = t120 * t169 + t121 * t165;
t309 = -t164 * t24 + t168 * t23;
t308 = -m(3) - m(5) - m(4);
t307 = qJD(3) + qJD(4);
t155 = pkin(5) * t235;
t140 = -qJD(2) * pkin(2) + t155;
t80 = -pkin(3) * t125 + t140;
t305 = -t80 * mrSges(5,1) + t18 * mrSges(5,3);
t304 = t80 * mrSges(5,2) - t17 * mrSges(5,3);
t202 = mrSges(3,1) * t169 - mrSges(3,2) * t165;
t303 = t165 * t330 + mrSges(2,1) + t202;
t296 = Ifges(4,1) * t290 + Ifges(4,4) * t289 + Ifges(4,5) * t283;
t277 = Ifges(5,4) * t64;
t26 = Ifges(5,2) * t205 + t144 * Ifges(5,6) + t277;
t295 = -t26 / 0.2e1;
t294 = t26 / 0.2e1;
t58 = Ifges(5,4) * t205;
t27 = Ifges(5,1) * t64 + Ifges(5,5) * t144 + t58;
t293 = -t27 / 0.2e1;
t292 = t27 / 0.2e1;
t260 = Ifges(4,4) * t126;
t51 = t125 * Ifges(4,2) + t148 * Ifges(4,6) + t260;
t291 = -t51 / 0.2e1;
t288 = -t205 / 0.2e1;
t286 = -t64 / 0.2e1;
t280 = -t144 / 0.2e1;
t276 = pkin(3) * t126;
t275 = pkin(3) * t164;
t271 = g(3) * t165;
t170 = cos(qJ(1));
t166 = sin(qJ(1));
t240 = t166 * t169;
t92 = t157 * t240 + t158 * t170;
t93 = t157 * t170 - t158 * t240;
t266 = -t92 * mrSges(5,1) + t93 * mrSges(5,2);
t238 = t169 * t170;
t94 = -t157 * t238 + t158 * t166;
t95 = t157 * t166 + t158 * t238;
t265 = t94 * mrSges(5,1) - t95 * mrSges(5,2);
t264 = mrSges(4,3) * t125;
t263 = mrSges(4,3) * t126;
t261 = Ifges(3,4) * t169;
t259 = Ifges(4,4) * t164;
t258 = Ifges(4,4) * t168;
t245 = t164 * t165;
t244 = t164 * t166;
t242 = t164 * t170;
t132 = t203 * qJD(2);
t233 = qJD(2) * t165;
t221 = pkin(5) * t233;
t237 = t168 * t132 + t164 * t221;
t229 = qJD(3) * t165;
t224 = Ifges(5,5) * t15 + Ifges(5,6) * t16 + Ifges(5,3) * t114;
t223 = Ifges(4,5) * t56 + Ifges(4,6) * t57 + Ifges(4,3) * t122;
t214 = t164 * t231;
t100 = -qJDD(2) * pkin(2) + t121;
t199 = mrSges(4,1) * t164 + mrSges(4,2) * t168;
t198 = -mrSges(5,1) * t157 - mrSges(5,2) * t158;
t197 = Ifges(4,1) * t168 - t259;
t196 = Ifges(4,1) * t164 + t258;
t194 = -Ifges(4,2) * t164 + t258;
t193 = Ifges(4,2) * t168 + t259;
t192 = Ifges(3,5) * t169 - Ifges(3,6) * t165;
t191 = Ifges(4,5) * t168 - Ifges(4,6) * t164;
t190 = Ifges(4,5) * t164 + Ifges(4,6) * t168;
t124 = t168 * t136;
t66 = -pkin(7) * t241 + t124 + (-pkin(5) * t164 - pkin(3)) * t169;
t73 = -pkin(7) * t245 + t83;
t33 = -t163 * t73 + t167 * t66;
t34 = t163 * t66 + t167 * t73;
t189 = t169 * t151 - t165 * t171;
t128 = t163 * t168 + t164 * t167;
t186 = t224 - t322;
t107 = -t164 * t238 + t166 * t168;
t105 = t164 * t240 + t168 * t170;
t182 = t128 * t169;
t181 = t188 * t169;
t178 = -t164 * t229 + t168 * t231;
t177 = t165 * t228 + t214;
t176 = Ifges(4,5) * t165 + t169 * t197;
t175 = Ifges(4,6) * t165 + t169 * t194;
t174 = Ifges(4,3) * t165 + t169 * t191;
t44 = t164 * t132 + t136 * t228 + (-t165 * t232 - t169 * t230) * pkin(5);
t71 = t307 * t128;
t138 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t234;
t135 = (pkin(5) + t275) * t165;
t116 = t199 * t165;
t108 = t168 * t238 + t244;
t106 = -t166 * t239 + t242;
t97 = t128 * t165;
t89 = qJD(1) * t181;
t88 = qJD(1) * t182;
t82 = -pkin(5) * t243 + t124;
t81 = pkin(3) * t177 + t220;
t79 = mrSges(4,1) * t148 - t263;
t78 = -mrSges(4,2) * t148 + t264;
t77 = -pkin(5) * t215 + t109;
t49 = mrSges(5,1) * t144 - mrSges(5,3) * t64;
t48 = -mrSges(5,2) * t144 + mrSges(5,3) * t205;
t45 = -qJD(3) * t83 + t237;
t43 = -mrSges(4,2) * t122 + mrSges(4,3) * t57;
t42 = mrSges(4,1) * t122 - mrSges(4,3) * t56;
t41 = -pkin(3) * t57 + t100;
t37 = -qJD(2) * t182 + t307 * t98;
t36 = -qJD(2) * t181 - t165 * t71;
t35 = -pkin(7) * t177 + t44;
t32 = t187 * qJD(2) + (-t149 + (pkin(7) * t165 - t136) * t164) * qJD(3) + t237;
t31 = -mrSges(5,1) * t205 + mrSges(5,2) * t64;
t30 = -mrSges(4,1) * t57 + mrSges(4,2) * t56;
t21 = t56 * Ifges(4,4) + t57 * Ifges(4,2) + t122 * Ifges(4,6);
t20 = t167 * t46 - t256;
t19 = -t163 * t46 - t252;
t11 = -mrSges(5,2) * t114 + mrSges(5,3) * t16;
t10 = mrSges(5,1) * t114 - mrSges(5,3) * t15;
t8 = -qJD(4) * t34 - t163 * t35 + t167 * t32;
t7 = qJD(4) * t33 + t163 * t32 + t167 * t35;
t6 = -mrSges(5,1) * t16 + mrSges(5,2) * t15;
t1 = [(m(4) * pkin(5) * t100 + Ifges(3,1) * t134 + Ifges(3,5) * qJDD(2) + t191 * t283 + t194 * t289 + t197 * t290) * t165 + t41 * (mrSges(5,1) * t97 - mrSges(5,2) * t98) - (t164 * t52 + t168 * t51) * t229 / 0.2e1 + t310 * mrSges(3,3) + (-qJDD(2) * mrSges(3,1) + mrSges(3,3) * t134 + t30) * pkin(5) * t165 + (-Ifges(5,1) * t98 - Ifges(5,4) * t97) * t298 + (Ifges(5,1) * t36 + Ifges(5,4) * t37) * t285 + (-t17 * t36 + t18 * t37 - t2 * t97 + t3 * t98) * mrSges(5,3) + (-t177 * t68 - t178 * t67 - t23 * t245 - t24 * t241) * mrSges(4,3) - pkin(1) * (-mrSges(3,1) * t133 + mrSges(3,2) * t134) + t135 * t6 + (-Ifges(5,5) * t98 - Ifges(5,6) * t97) * t284 + (Ifges(5,5) * t36 + Ifges(5,6) * t37) * t279 + t100 * t116 + t44 * t78 + t45 * t79 + t80 * (-mrSges(5,1) * t37 + mrSges(5,2) * t36) + t81 * t31 + t82 * t42 + t83 * t43 + t7 * t48 + t8 * t49 + t33 * t10 + t34 * t11 + t140 * (mrSges(4,1) * t177 + mrSges(4,2) * t178) + (t262 + t195) * t133 / 0.2e1 + t261 * t328 - t138 * t221 + t202 * t249 + (Ifges(3,4) * t328 - t224 / 0.2e1 - t223 / 0.2e1 - Ifges(5,3) * t284 - Ifges(5,6) * t297 - Ifges(5,5) * t298 + (-Ifges(3,2) * t165 + t261) * t227 / 0.2e1 - Ifges(4,6) * t289 - Ifges(4,5) * t290 - Ifges(4,3) * t283 + (-mrSges(3,2) * pkin(5) + Ifges(3,6)) * qJDD(2) + (pkin(5) * mrSges(3,3) + t329) * t133 + t321 + t322) * t169 + m(4) * (t140 * t220 + t23 * t83 + t24 * t82 + t44 * t68 + t45 * t67) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t310) + t312 * t231 / 0.2e1 - t98 * t300 - t97 * t301 + t214 * t291 + t36 * t292 + t37 * t294 + t241 * t296 + (qJD(2) * t176 - t196 * t229) * t281 + (-Ifges(5,4) * t98 - Ifges(5,2) * t97) * t297 + (Ifges(5,4) * t36 + Ifges(5,2) * t37) * t287 + (qJD(2) * t174 - t190 * t229) * t332 + (qJD(2) * t175 - t193 * t229) * t333 + qJD(2) ^ 2 * t192 / 0.2e1 + (t320 + t334) * t233 + t323 * t227 - t313 * t220 - t21 * t245 / 0.2e1 + (-t244 * t299 - t108 * mrSges(4,1) - t95 * mrSges(5,1) - t107 * mrSges(4,2) - t94 * mrSges(5,2) + t308 * (t170 * pkin(1) + t166 * pkin(5)) + t318 * t166 + (-m(4) * t204 - m(5) * t189 - t303) * t170) * g(2) + Ifges(2,3) * qJDD(1) + (-t242 * t299 - t106 * mrSges(4,1) - t93 * mrSges(5,1) - t105 * mrSges(4,2) - t92 * mrSges(5,2) + (m(3) * pkin(1) - m(5) * (-pkin(1) - t189) - m(4) * t136 + t303) * t166 + (t308 * pkin(5) + t318) * t170) * g(1) + m(5) * (t135 * t41 + t17 * t8 + t18 * t7 + t2 * t34 + t3 * t33 + t80 * t81); (t125 * t194 + t126 * t197 + t148 * t191) * qJD(3) / 0.2e1 + (-t17 * t89 + t18 * t88) * mrSges(5,3) - (Ifges(5,4) * t285 + Ifges(5,2) * t287 + Ifges(5,6) * t279 + t294 + t305) * t71 + t168 * t21 / 0.2e1 - t151 * t6 + Ifges(3,6) * t133 + Ifges(3,5) * t134 - t120 * mrSges(3,2) - t121 * mrSges(3,1) - (t154 + t312) * t234 / 0.2e1 - t77 * t78 - t76 * t79 + t74 * t10 + t75 * t11 - pkin(2) * t30 - t192 * t227 / 0.2e1 + (-Ifges(5,1) * t89 - Ifges(5,4) * t88) * t286 + (-Ifges(5,5) * t89 - Ifges(5,6) * t88) * t280 - t80 * (mrSges(5,1) * t88 - mrSges(5,2) * t89) + (-Ifges(5,4) * t89 - Ifges(5,2) * t88) * t288 + t52 * t228 / 0.2e1 + (-pkin(2) * t100 - t140 * t156 - t67 * t76 - t68 * t77) * m(4) + t51 * t216 / 0.2e1 + (-(Ifges(5,1) * t285 + Ifges(5,4) * t287 + Ifges(5,5) * t279 + t292 + t304) * t307 - t2 * mrSges(5,3) + t41 * mrSges(5,1) - 0.2e1 * t301) * t188 + (t41 * mrSges(5,2) - mrSges(5,3) * t3 + 0.2e1 * t300) * t128 + t138 * t155 + (Ifges(5,5) * t286 + Ifges(5,6) * t288 + Ifges(5,3) * t280 + t234 * t329 - t320) * t235 + (-t228 * t67 - t230 * t68 + (t239 * t67 + t243 * t68) * qJD(1) + t309) * mrSges(4,3) + (t168 * t43 - t164 * t42 - t79 * t228 - t78 * t230 + m(4) * ((-t164 * t68 - t168 * t67) * qJD(3) + t309)) * pkin(6) + t311 * t31 + (t125 * t175 + t126 * t176 + t148 * t174) * t331 + t190 * t283 + t193 * t289 + t196 * t290 + t230 * t291 - t89 * t293 - t88 * t295 + t164 * t296 + t100 * t200 + (g(1) * t170 + g(2) * t166) * (t165 * t325 + t169 * t324 + t201) + (t165 * t324 - t169 * t325 - t202) * g(3) - t323 * qJD(1) ^ 2 + t148 * t140 * t199 + t313 * t156 + t316 * t49 + (-t151 * t41 + t17 * t316 + t18 * t317 + t2 * t75 + t3 * t74 + t311 * t80) * m(5) + t317 * t48 + Ifges(3,3) * qJDD(2); -t321 - (-Ifges(4,2) * t126 + t119 + t52) * t125 / 0.2e1 + (Ifges(5,1) * t286 + Ifges(5,4) * t288 + Ifges(5,5) * t280 + t293 - t304) * t205 - (Ifges(5,4) * t286 + Ifges(5,2) * t288 + Ifges(5,6) * t280 + t295 - t305) * t64 - t148 * (Ifges(4,5) * t125 - Ifges(4,6) * t126) / 0.2e1 - t140 * (mrSges(4,1) * t126 + mrSges(4,2) * t125) + (t79 + t263) * t68 + t223 + g(3) * t116 - t20 * t48 - t19 * t49 + (-t78 + t264) * t67 + (t163 * t2 + t167 * t3 + (-t163 * t17 + t167 * t18) * qJD(4)) * t299 + t51 * t281 + t186 + (-mrSges(4,2) * t106 + t105 * t314 - t266) * g(2) + (mrSges(4,2) * t108 - t107 * t314 - t265) * g(1) - t126 * (Ifges(4,1) * t125 - t260) / 0.2e1 - (-m(5) * t275 + t198) * t271 - t31 * t276 - m(5) * (t17 * t19 + t18 * t20 + t276 * t80) + ((-t163 * t49 + t167 * t48) * qJD(4) + t10 * t167 + t11 * t163) * pkin(3); -t80 * (mrSges(5,1) * t64 + mrSges(5,2) * t205) + (Ifges(5,1) * t205 - t277) * t286 + t26 * t285 + (Ifges(5,5) * t205 - Ifges(5,6) * t64) * t280 - t17 * t48 + t18 * t49 - g(1) * t265 - g(2) * t266 - t198 * t271 + (t17 * t205 + t18 * t64) * mrSges(5,3) + t186 + (-Ifges(5,2) * t64 + t27 + t58) * t288;];
tau = t1;

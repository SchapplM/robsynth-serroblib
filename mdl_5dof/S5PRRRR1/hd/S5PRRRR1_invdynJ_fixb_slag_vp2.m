% Calculate vector of inverse dynamics joint torques for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:21
% EndTime: 2019-07-18 13:28:30
% DurationCPUTime: 3.92s
% Computational Cost: add. (2539->395), mult. (6328->537), div. (0->0), fcn. (4698->10), ass. (0->194)
t130 = qJD(3) + qJD(4);
t279 = t130 * Ifges(5,6) / 0.2e1;
t140 = cos(qJ(3));
t201 = qJD(2) * t140;
t141 = cos(qJ(2));
t204 = qJD(1) * t141;
t113 = -pkin(2) * t201 - t204;
t223 = t130 * Ifges(5,5);
t278 = t113 * mrSges(5,2) + t223 / 0.2e1;
t139 = cos(qJ(4));
t135 = sin(qJ(4));
t136 = sin(qJ(3));
t211 = t135 * t136;
t160 = -t139 * t140 + t211;
t101 = t160 * qJD(2);
t210 = t135 * t140;
t161 = t136 * t139 + t210;
t102 = t161 * qJD(2);
t224 = t102 * Ifges(5,4);
t277 = t279 + t224 / 0.2e1 - t101 * Ifges(5,2) / 0.2e1;
t137 = sin(qJ(2));
t274 = t130 * t137;
t133 = qJ(3) + qJ(4);
t127 = sin(t133);
t128 = cos(t133);
t134 = sin(qJ(5));
t235 = mrSges(6,2) * t134;
t273 = mrSges(6,3) * t128 + t127 * t235;
t138 = cos(qJ(5));
t197 = qJD(1) * qJD(2);
t186 = t141 * t197;
t110 = qJDD(1) * t137 + t186;
t205 = qJD(1) * t137;
t185 = qJD(3) * t205;
t84 = -t110 * t136 - t140 * t185;
t147 = qJDD(3) * pkin(2) + t84;
t157 = qJD(3) * pkin(2) - t136 * t205;
t188 = t140 * t205;
t81 = t135 * t188 - t139 * t157;
t218 = qJD(4) * t81;
t83 = t110 * t140 - t136 * t185;
t25 = t135 * t147 + t139 * t83 - t218;
t82 = t135 * t157 + t139 * t188;
t53 = t113 * t138 - t134 * t82;
t196 = qJD(2) * qJD(3);
t107 = qJDD(2) * t140 - t136 * t196;
t109 = -t141 * qJDD(1) + t137 * t197;
t80 = -pkin(2) * t107 + t109;
t7 = qJD(5) * t53 + t134 * t80 + t138 * t25;
t54 = t113 * t134 + t138 * t82;
t8 = -qJD(5) * t54 - t134 * t25 + t138 * t80;
t272 = -t8 * t134 + t138 * t7;
t265 = t161 * qJD(4);
t271 = qJD(3) * t161 + t265;
t108 = qJDD(2) * t136 + t140 * t196;
t50 = -qJD(2) * t265 + t107 * t139 - t108 * t135;
t166 = t134 * t54 + t138 * t53;
t170 = Ifges(6,5) * t138 - Ifges(6,6) * t134;
t229 = Ifges(6,4) * t138;
t172 = -Ifges(6,2) * t134 + t229;
t230 = Ifges(6,4) * t134;
t174 = Ifges(6,1) * t138 - t230;
t251 = t138 / 0.2e1;
t252 = -t134 / 0.2e1;
t79 = t102 * t138 + t130 * t134;
t257 = t79 / 0.2e1;
t244 = t79 * Ifges(6,4);
t78 = -t102 * t134 + t130 * t138;
t98 = qJD(5) + t101;
t33 = t78 * Ifges(6,2) + t98 * Ifges(6,6) + t244;
t75 = Ifges(6,4) * t78;
t34 = t79 * Ifges(6,1) + t98 * Ifges(6,5) + t75;
t270 = -t166 * mrSges(6,3) + t98 * t170 / 0.2e1 + t174 * t257 + t78 * t172 / 0.2e1 + t34 * t251 + t33 * t252;
t269 = -t113 * mrSges(5,1) - t53 * mrSges(6,1) + t54 * mrSges(6,2) + t277;
t268 = -t33 / 0.2e1;
t267 = -t134 * t53 + t138 * t54;
t129 = qJDD(3) + qJDD(4);
t153 = t160 * qJD(4);
t49 = -qJD(2) * t153 + t107 * t135 + t108 * t139;
t23 = qJD(5) * t78 + t129 * t134 + t138 * t49;
t24 = -qJD(5) * t79 + t129 * t138 - t134 * t49;
t48 = qJDD(5) - t50;
t264 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t48;
t263 = -mrSges(5,1) * t128 + t127 * (mrSges(5,2) - mrSges(6,3)) - mrSges(3,1);
t262 = t23 / 0.2e1;
t261 = t24 / 0.2e1;
t260 = t48 / 0.2e1;
t259 = -t78 / 0.2e1;
t258 = -t79 / 0.2e1;
t256 = -t98 / 0.2e1;
t255 = m(5) + m(6);
t254 = t101 / 0.2e1;
t253 = -t102 / 0.2e1;
t250 = t140 / 0.2e1;
t249 = -mrSges(5,1) * t129 - mrSges(6,1) * t24 + mrSges(6,2) * t23 + mrSges(5,3) * t49;
t248 = g(1) * t141;
t247 = g(2) * t127;
t245 = t78 * Ifges(6,6);
t243 = t79 * Ifges(6,5);
t87 = t160 * t205;
t241 = t81 * t87;
t88 = t161 * t204;
t240 = t81 * t88;
t239 = t98 * Ifges(6,3);
t233 = mrSges(5,3) * t102;
t237 = -mrSges(5,1) * t130 - mrSges(6,1) * t78 + mrSges(6,2) * t79 + t233;
t236 = mrSges(6,1) * t138;
t234 = mrSges(5,3) * t101;
t231 = Ifges(4,4) * t136;
t97 = Ifges(5,4) * t101;
t228 = Ifges(4,5) * t140;
t227 = Ifges(4,6) * t136;
t225 = t102 * Ifges(5,1);
t219 = t140 * Ifges(4,2);
t217 = qJD(4) * t82;
t216 = t101 * t134;
t215 = t101 * t138;
t214 = t109 * t141;
t213 = t134 * t137;
t212 = t134 * t141;
t209 = t137 * t138;
t208 = t138 * t141;
t207 = t273 * t137;
t206 = t273 * t141;
t203 = qJD(2) * t136;
t202 = qJD(2) * t137;
t200 = qJD(2) * t141;
t199 = qJD(5) * t134;
t198 = qJD(5) * t138;
t195 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t194 = t136 ^ 2 + t140 ^ 2 - 0.1e1;
t192 = pkin(2) * t203;
t66 = mrSges(5,1) * t101 + mrSges(5,2) * t102;
t187 = m(5) * t113 + t66;
t183 = mrSges(5,1) + t236;
t181 = g(3) * t137 + t248;
t180 = g(1) * t137 - g(3) * t141;
t179 = -t7 * t134 - t8 * t138;
t26 = t135 * t83 - t139 * t147 + t217;
t37 = t200 * t210 + (t136 * t200 + t140 * t274) * t139 - t211 * t274;
t91 = t161 * t137;
t178 = t26 * t91 + t37 * t81;
t118 = -mrSges(4,1) * t140 + mrSges(4,2) * t136;
t177 = mrSges(4,1) * t136 + mrSges(4,2) * t140;
t176 = -t235 + t236;
t175 = mrSges(6,1) * t134 + mrSges(6,2) * t138;
t173 = Ifges(6,1) * t134 + t229;
t171 = Ifges(6,2) * t138 + t230;
t169 = Ifges(6,5) * t134 + Ifges(6,6) * t138;
t51 = -mrSges(6,2) * t98 + mrSges(6,3) * t78;
t52 = mrSges(6,1) * t98 - mrSges(6,3) * t79;
t168 = t134 * t52 - t138 * t51;
t167 = -t134 * t51 - t138 * t52;
t92 = t160 * t137;
t76 = t134 * t92 - t208;
t164 = t138 * t92 + t212;
t116 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t203;
t117 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t201;
t162 = t116 * t140 + t117 * t136;
t159 = mrSges(5,3) + t175;
t85 = -mrSges(5,2) * t130 - t234;
t158 = -t168 + t85;
t156 = t136 * (Ifges(4,1) * t140 - t231);
t149 = -t127 * mrSges(5,2) + (t183 - t235) * t128;
t148 = mrSges(5,2) * t128 + t127 * t183;
t146 = t148 + t177;
t124 = Ifges(4,4) * t201;
t145 = -t136 * (Ifges(4,6) * qJD(3) + (t219 + t231) * qJD(2)) / 0.2e1 + (Ifges(4,1) * t203 + Ifges(4,5) * qJD(3) + t124) * t250 - t177 * t204;
t3 = t23 * Ifges(6,4) + t24 * Ifges(6,2) + t48 * Ifges(6,6);
t32 = t239 + t243 + t245;
t4 = t23 * Ifges(6,1) + t24 * Ifges(6,4) + t48 * Ifges(6,5);
t60 = t223 - t97 + t225;
t144 = t216 * t268 + t134 * t4 / 0.2e1 + Ifges(5,3) * t129 + t169 * t260 + t171 * t261 + t173 * t262 + t3 * t251 + t81 * t234 + Ifges(5,5) * t49 + Ifges(5,6) * t50 - t25 * mrSges(5,2) + t34 * t215 / 0.2e1 + (-mrSges(5,1) - t176) * t26 + (-t97 + t60) * t254 + (t32 - t224) * t253 + (-Ifges(5,1) * t253 - t170 * t256 - t172 * t259 - t174 * t258 + t278) * t101 + (Ifges(6,5) * t258 - Ifges(5,2) * t254 + Ifges(6,6) * t259 + Ifges(6,3) * t256 + t269 + t279) * t102 + (-t53 * t215 - t54 * t216 + t247 + t272) * mrSges(6,3) + (t81 * t175 + t270) * qJD(5);
t142 = qJD(2) ^ 2;
t106 = t118 * qJD(2);
t96 = t128 * t208 + t213;
t95 = -t128 * t212 + t209;
t94 = -t128 * t209 + t212;
t93 = t128 * t213 + t208;
t90 = t160 * t204;
t89 = t161 * t205;
t70 = t134 * t205 - t138 * t90;
t69 = t134 * t90 + t138 * t205;
t68 = t134 * t192 - t138 * t89;
t67 = t134 * t89 + t138 * t192;
t61 = t175 * t101;
t42 = -mrSges(5,2) * t129 + mrSges(5,3) * t50;
t36 = -t141 * t101 - t271 * t137;
t13 = -mrSges(5,1) * t50 + mrSges(5,2) * t49;
t12 = qJD(5) * t164 - t134 * t36 + t138 * t202;
t11 = qJD(5) * t76 + t134 * t202 + t138 * t36;
t10 = -mrSges(6,2) * t48 + mrSges(6,3) * t24;
t9 = mrSges(6,1) * t48 - mrSges(6,3) * t23;
t1 = [m(2) * qJDD(1) - t164 * t10 + t11 * t51 + t12 * t52 + t36 * t85 - t92 * t42 + t76 * t9 + t249 * t91 + t237 * t37 + (-m(2) - m(3) - m(4) - t255) * g(3) + (qJDD(2) * mrSges(3,1) + t107 * mrSges(4,1) - t142 * mrSges(3,2) - t108 * mrSges(4,2) - t13 + (-t116 * t136 + t117 * t140) * qJD(2)) * t141 + (t140 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t107) - t142 * mrSges(3,1) - t136 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t108) - qJDD(2) * mrSges(3,2) - t162 * qJD(3) + (t106 + t66) * qJD(2)) * t137 + m(4) * (-t214 + (-t136 * t84 + t140 * t83 + t186 * t194) * t137) + m(3) * (t110 * t137 - t214) + m(5) * (t113 * t202 - t141 * t80 - t25 * t92 + t36 * t82 + t178) + m(6) * (t11 * t54 + t12 * t53 - t164 * t7 + t76 * t8 + t178); -t110 * mrSges(3,2) + t90 * t85 - t69 * t52 - t70 * t51 - g(3) * (mrSges(6,1) * t96 + mrSges(6,2) * t95) - g(1) * (mrSges(6,1) * t94 + mrSges(6,2) * t93) + Ifges(3,3) * qJDD(2) - t237 * t88 + (-mrSges(3,1) + t118) * t109 - m(6) * (t53 * t69 + t54 * t70 + t240) - m(5) * (-t82 * t90 + t240) + (t239 / 0.2e1 + t245 / 0.2e1 + t243 / 0.2e1 + t32 / 0.2e1 - t82 * mrSges(5,3) - t269 - t277) * t271 + (-t97 / 0.2e1 + t225 / 0.2e1 + t60 / 0.2e1 + t159 * t81 + t270 + t278) * (-qJD(3) * t160 - t153) + (t195 * g(3) - t263 * g(1) + (qJD(2) * mrSges(3,1) - t106 - t187) * qJD(1)) * t137 + (-m(4) * t194 * qJD(1) ^ 2 * t137 + mrSges(3,2) * t197 + t195 * g(1) + t263 * g(3)) * t141 + (t180 * mrSges(4,1) + t83 * mrSges(4,3) + Ifges(4,4) * t108 + Ifges(4,2) * t107 + Ifges(4,6) * qJDD(3) - t117 * t204) * t140 + (-t180 * mrSges(4,2) - t84 * mrSges(4,3) + Ifges(4,1) * t108 + Ifges(4,4) * t107 + Ifges(4,5) * qJDD(3) + t116 * t204) * t136 + ((t228 / 0.2e1 - t227 / 0.2e1) * qJD(3) + (t156 / 0.2e1 + (Ifges(4,4) * t140 - Ifges(4,2) * t136) * t250) * qJD(2) + t145) * qJD(3) + (t80 * mrSges(5,1) - t25 * mrSges(5,3) - Ifges(5,4) * t49 - Ifges(5,2) * t50 - Ifges(5,6) * t129 + t264) * t160 + (t170 * t260 + t174 * t262 + t172 * t261 + Ifges(5,5) * t129 + Ifges(5,1) * t49 + Ifges(5,4) * t50 + t4 * t251 + t3 * t252 + t80 * mrSges(5,2) + t159 * t26 + t179 * mrSges(6,3) + (-mrSges(6,3) * t267 + t138 * t268 + t169 * t256 + t171 * t259 + t173 * t258 + t81 * t176 + t34 * t252) * qJD(5)) * t161 + ((m(6) * t166 - t167 + t187) * t136 * qJD(3) + (-t134 * t10 - t138 * t9 - t13 + t168 * qJD(5) + (-t198 * t54 + t199 * t53 + t179 + t180) * m(6) + (t180 - t80) * m(5)) * t140) * pkin(2); (g(3) * t146 + qJD(1) * t162) * t137 + (-t118 + t149) * g(2) + t237 * t87 + t82 * t233 - m(5) * (-t82 * t89 - t241) + Ifges(4,5) * t108 + Ifges(4,6) * t107 + t89 * t85 + t81 * t61 - t83 * mrSges(4,2) + t84 * mrSges(4,1) - t67 * t52 - t68 * t51 + (-qJD(3) * (-t227 + t228) / 0.2e1 - t140 * t124 / 0.2e1 + (t136 * t219 / 0.2e1 - t156 / 0.2e1) * qJD(2) - t145) * qJD(2) - m(6) * (t53 * t67 + t54 * t68 - t241) + (t255 * t140 * g(2) + (-qJD(2) * t66 + t181 * m(6) + (-qJD(2) * t113 + t181) * m(5)) * t136 + (-m(6) * t26 + m(5) * (-t26 + t217) - t249 + (m(6) * t267 + t158) * qJD(4)) * t139 + (t138 * t10 - t134 * t9 + t42 + t167 * qJD(5) + t237 * qJD(4) + m(6) * (-t198 * t53 - t199 * t54 + t218 + t272) + m(5) * (t25 + t218)) * t135) * pkin(2) - g(1) * t206 - g(3) * t207 + t144 + Ifges(4,3) * qJDD(3) + t146 * t248; (t233 - t237) * t82 + (t141 * t148 - t206) * g(1) + t149 * g(2) + (t137 * t148 - t207) * g(3) + (-m(6) * (-t267 + t82) + t61 + t158) * t81 + t144; -t81 * (mrSges(6,1) * t79 + mrSges(6,2) * t78) + (Ifges(6,1) * t78 - t244) * t258 + t33 * t257 + (Ifges(6,5) * t78 - Ifges(6,6) * t79) * t256 - t53 * t51 + t54 * t52 - g(1) * (mrSges(6,1) * t95 - mrSges(6,2) * t96) - g(3) * (-mrSges(6,1) * t93 + mrSges(6,2) * t94) - t175 * t247 + (t53 * t78 + t54 * t79) * mrSges(6,3) + (-Ifges(6,2) * t79 + t34 + t75) * t259 + t264;];
tau  = t1;

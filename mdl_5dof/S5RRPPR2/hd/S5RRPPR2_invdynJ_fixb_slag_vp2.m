% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:04
% EndTime: 2019-12-05 18:20:10
% DurationCPUTime: 2.74s
% Computational Cost: add. (2409->278), mult. (3834->388), div. (0->0), fcn. (2190->14), ass. (0->153)
t133 = sin(qJ(5));
t254 = -t133 / 0.2e1;
t129 = sin(pkin(9));
t125 = t129 ^ 2;
t131 = cos(pkin(9));
t253 = mrSges(5,3) * (t131 ^ 2 + t125);
t245 = t131 * mrSges(5,1) - mrSges(5,2) * t129;
t252 = mrSges(4,1) + t245;
t127 = qJD(1) + qJD(2);
t130 = sin(pkin(8));
t132 = cos(pkin(8));
t134 = sin(qJ(2));
t214 = qJD(1) * pkin(1);
t185 = t134 * t214;
t137 = cos(qJ(2));
t184 = t137 * t214;
t95 = pkin(2) * t127 + t184;
t53 = t130 * t95 + t132 * t185;
t47 = qJ(4) * t127 + t53;
t37 = -t131 * qJD(3) + t129 * t47;
t209 = t129 * t37;
t38 = qJD(3) * t129 + t131 * t47;
t239 = t131 * t38;
t159 = t209 + t239;
t236 = t127 * t253;
t251 = -m(5) * t159 - t236;
t250 = mrSges(4,2) - mrSges(5,3);
t156 = -pkin(4) * t131 - pkin(7) * t129 - pkin(3);
t249 = -m(6) * t156 + t129 * mrSges(6,3) + t252;
t100 = -t127 * t131 + qJD(5);
t136 = cos(qJ(5));
t215 = Ifges(6,4) * t136;
t147 = (-t133 * Ifges(6,2) + t215) * t129;
t216 = Ifges(6,4) * t133;
t148 = (Ifges(6,1) * t136 - t216) * t129;
t101 = t130 * t185;
t52 = t132 * t95 - t101;
t166 = qJD(4) - t52;
t25 = t127 * t156 + t166;
t8 = -t133 * t38 + t136 * t25;
t9 = t133 * t25 + t136 * t38;
t167 = -t133 * t8 + t136 * t9;
t248 = t37 * (mrSges(6,1) * t136 - mrSges(6,2) * t133) - t167 * mrSges(6,3) + (-t136 * t147 / 0.2e1 + t148 * t254) * t127 + (0.2e1 * Ifges(6,5) * t254 - Ifges(6,6) * t136) * t100;
t124 = qJDD(1) + qJDD(2);
t192 = qJD(4) * t127;
t231 = pkin(1) * t134;
t183 = qJD(2) * t231;
t229 = pkin(1) * t137;
t92 = -qJD(1) * t183 + qJDD(1) * t229;
t75 = pkin(2) * t124 + t92;
t193 = qJD(2) * t137;
t93 = (qJD(1) * t193 + qJDD(1) * t134) * pkin(1);
t36 = t130 * t75 + t132 * t93;
t22 = qJ(4) * t124 + t192 + t36;
t16 = -t131 * qJDD(3) + t129 * t22;
t17 = qJDD(3) * t129 + t131 * t22;
t204 = t131 * t17;
t246 = t129 * t16 + t204;
t243 = (-Ifges(6,2) * t136 - t216) * t254 + t136 * (-Ifges(6,1) * t133 - t215) / 0.2e1;
t191 = qJD(4) * t131;
t196 = t131 * t133;
t225 = pkin(2) * t130;
t109 = qJ(4) + t225;
t195 = t131 * t136;
t224 = pkin(2) * t132;
t89 = t156 - t224;
t44 = t109 * t195 + t133 * t89;
t194 = t132 * t134;
t152 = pkin(1) * (t130 * t137 + t194);
t78 = qJD(1) * t152;
t80 = t132 * t184 - t101;
t241 = -qJD(5) * t44 - t133 * t191 - t136 * t78 + t196 * t80;
t43 = -t109 * t196 + t136 * t89;
t240 = qJD(5) * t43 - t133 * t78 + t136 * t191 - t195 * t80;
t237 = (mrSges(3,1) * t231 + mrSges(3,2) * t229) * t127;
t190 = qJD(5) * t127;
t70 = (t124 * t136 - t133 * t190) * t129;
t71 = (-t124 * t133 - t136 * t190) * t129;
t23 = -mrSges(6,1) * t71 + mrSges(6,2) * t70;
t235 = t124 * t253 + t129 * t23;
t128 = qJ(1) + qJ(2);
t120 = pkin(8) + t128;
t111 = sin(t120);
t112 = cos(t120);
t121 = sin(t128);
t122 = cos(t128);
t66 = -t111 * t136 + t112 * t196;
t67 = -t111 * t133 - t112 * t195;
t234 = t122 * mrSges(3,1) - t67 * mrSges(6,1) - t121 * mrSges(3,2) - t66 * mrSges(6,2) - t250 * t111 + t249 * t112;
t64 = t111 * t196 + t112 * t136;
t65 = t111 * t195 - t112 * t133;
t233 = t121 * mrSges(3,1) + t65 * mrSges(6,1) + t122 * mrSges(3,2) - t64 * mrSges(6,2) + t249 * t111 + t250 * t112;
t232 = m(4) * t52 - m(5) * t166 + (m(5) * pkin(3) + t252) * t127;
t135 = sin(qJ(1));
t230 = pkin(1) * t135;
t138 = cos(qJ(1));
t228 = pkin(1) * t138;
t227 = pkin(2) * t121;
t226 = pkin(2) * t122;
t221 = mrSges(4,1) * t124;
t220 = mrSges(4,2) * t124;
t212 = t127 * mrSges(4,2);
t161 = mrSges(6,1) * t133 + mrSges(6,2) * t136;
t199 = t127 * t129;
t69 = t161 * t199;
t208 = t129 * t69;
t207 = t129 * t80;
t205 = t131 * t16;
t114 = pkin(2) + t229;
t84 = pkin(1) * t194 + t130 * t114;
t201 = t124 * t129;
t200 = t124 * t131;
t198 = t129 * t133;
t197 = t129 * t136;
t98 = qJDD(5) - t200;
t188 = Ifges(6,5) * t70 + Ifges(6,6) * t71 + Ifges(6,3) * t98;
t187 = mrSges(6,3) * t198;
t186 = mrSges(6,3) * t197;
t104 = t112 * qJ(4);
t174 = t104 - t227;
t35 = -t130 * t93 + t132 * t75;
t83 = t114 * t132 - t130 * t231;
t82 = -mrSges(5,1) * t200 + mrSges(5,2) * t201;
t169 = -t227 - t230;
t168 = -g(2) * t112 - g(3) * t111;
t165 = qJDD(4) - t35;
t160 = -t111 * qJ(4) - t226;
t62 = -mrSges(6,2) * t100 - t127 * t187;
t63 = mrSges(6,1) * t100 - t127 * t186;
t158 = t133 * t63 - t136 * t62;
t81 = t132 * pkin(1) * t193 - t130 * t183;
t155 = -pkin(3) * t111 + t174;
t51 = t156 - t83;
t76 = qJ(4) + t84;
t21 = t133 * t51 + t195 * t76;
t20 = t136 * t51 - t196 * t76;
t72 = qJD(4) + t81;
t154 = (t16 * t76 + t37 * t72) * t129;
t145 = -t112 * pkin(3) + t160;
t24 = -pkin(3) * t124 + t165;
t15 = t124 * t156 + t165;
t3 = qJD(5) * t8 + t133 * t15 + t136 * t17;
t4 = -qJD(5) * t9 - t133 * t17 + t136 * t15;
t88 = t161 * t129;
t139 = t98 * (-Ifges(6,3) * t131 + (Ifges(6,5) * t136 - Ifges(6,6) * t133) * t129) / 0.2e1 + (Ifges(5,4) * t129 + Ifges(5,2) * t131) * t200 + (Ifges(5,1) * t129 + Ifges(5,4) * t131) * t201 + (Ifges(6,1) * t70 + Ifges(6,4) * t71 + Ifges(6,5) * t98) * t197 / 0.2e1 - (Ifges(6,4) * t70 + Ifges(6,2) * t71 + Ifges(6,6) * t98) * t198 / 0.2e1 - t131 * t188 / 0.2e1 + t4 * (-mrSges(6,1) * t131 - t186) + t3 * (mrSges(6,2) * t131 - t187) - t24 * t245 + t71 * (-Ifges(6,6) * t131 + t147) / 0.2e1 + t70 * (-Ifges(6,5) * t131 + t148) / 0.2e1 + t16 * t88 + t92 * mrSges(3,1) - t93 * mrSges(3,2) + t35 * mrSges(4,1) - t36 * mrSges(4,2) + t243 * t125 * t190 + (Ifges(3,3) + Ifges(4,3)) * t124 + t246 * mrSges(5,3) + t248 * qJD(5) * t129;
t113 = -pkin(3) - t224;
t79 = qJD(2) * t152;
t77 = -pkin(3) - t83;
t40 = -mrSges(6,2) * t98 + mrSges(6,3) * t71;
t39 = mrSges(6,1) * t98 - mrSges(6,3) * t70;
t6 = -qJD(5) * t21 + t136 * t79 - t196 * t72;
t5 = qJD(5) * t20 + t133 * t79 + t195 * t72;
t1 = [-t84 * t220 - t81 * t212 + m(6) * (t20 * t4 + t21 * t3 + t5 * t9 + t6 * t8 + t154) + m(5) * (t24 * t77 + t154) + m(4) * (t35 * t83 + t36 * t84 + t53 * t81) + t139 + t77 * t82 + t5 * t62 + t6 * t63 + t20 * t39 + t21 * t40 + t83 * t221 + m(3) * (t134 * t93 + t137 * t92) * pkin(1) + Ifges(2,3) * qJDD(1) - t232 * t79 + (m(5) * t204 + t235) * t76 + (m(5) * t239 + t208 + t236) * t72 + (mrSges(3,1) * t229 - mrSges(3,2) * t231) * t124 - t237 * qJD(2) + (m(3) * t230 - m(5) * (t155 - t230) - m(6) * (t104 + t169) - m(4) * t169 + t135 * mrSges(2,1) + t138 * mrSges(2,2) + t233) * g(3) + (m(3) * t228 - m(4) * (-t226 - t228) - m(6) * (t160 - t228) - m(5) * (t145 - t228) + t138 * mrSges(2,1) - t135 * mrSges(2,2) + t234) * g(2); -t220 * t225 - t69 * t207 + m(5) * (t159 * qJD(4) + t113 * t24) + t139 + t113 * t82 + t43 * t39 + t44 * t40 + qJD(4) * t208 + t221 * t224 + m(4) * (t130 * t36 + t132 * t35) * pkin(2) + t241 * t63 + t240 * t62 + t192 * t253 + t237 * qJD(1) + (m(5) * t246 + t235) * t109 + (t3 * t44 + t4 * t43 + (qJD(4) * t37 + t109 * t16) * t129 - t207 * t37 + t240 * t9 + t241 * t8) * m(6) + t232 * t78 + (-m(4) * t53 + t212 + t251) * t80 + (m(4) * t227 - m(5) * t155 - m(6) * t174 + t233) * g(3) + (m(4) * t226 - m(5) * t145 - m(6) * t160 + t234) * g(2); m(4) * qJDD(3) - t131 * t23 + (-t133 * t39 + t136 * t40 + (-t133 * t62 - t136 * t63) * qJD(5)) * t129 + m(5) * (t129 * t17 - t205) + m(6) * (-t205 + (-t133 * t4 + t136 * t3 + (-t133 * t9 - t136 * t8) * qJD(5)) * t129) + (-m(4) - m(5) - m(6)) * g(1); t133 * t40 + t136 * t39 - t158 * qJD(5) + (t167 * qJD(5) + t133 * t3 + t136 * t4 + t168) * m(6) + (t168 + t24) * m(5) + t82 + (-t208 + t158 * t131 - m(6) * (t195 * t9 - t196 * t8 + t209) + t251) * t127; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t8 * t62 + t9 * t63 + g(1) * t88 - g(2) * (mrSges(6,1) * t64 + mrSges(6,2) * t65) - g(3) * (-mrSges(6,1) * t66 + mrSges(6,2) * t67) + (-t243 * t199 - t248) * t199 + t188;];
tau = t1;

% Calculate vector of inverse dynamics joint torques for
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:38
% EndTime: 2019-12-05 15:46:49
% DurationCPUTime: 4.55s
% Computational Cost: add. (2176->324), mult. (4520->445), div. (0->0), fcn. (3085->14), ass. (0->167)
t224 = m(5) + m(6) + m(4);
t226 = pkin(2) * t224;
t123 = sin(qJ(4));
t126 = cos(qJ(4));
t150 = -mrSges(5,1) * t126 + mrSges(5,2) * t123;
t225 = -mrSges(4,1) + t150;
t102 = pkin(4) * t126 + pkin(3);
t117 = qJ(4) + qJ(5);
t111 = sin(t117);
t112 = cos(t117);
t223 = m(5) * pkin(3) + m(6) * t102 + t112 * mrSges(6,1) - t111 * mrSges(6,2) - t225;
t222 = -m(5) * pkin(6) + m(6) * (-pkin(7) - pkin(6)) + mrSges(4,2) - mrSges(6,3);
t164 = qJD(2) * qJD(4);
t89 = qJDD(2) * t126 - t123 * t164;
t221 = t89 / 0.2e1;
t116 = qJ(2) + pkin(9);
t105 = sin(t116);
t201 = g(3) * t105;
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t118 = sin(pkin(9));
t206 = pkin(2) * t118;
t99 = pkin(6) + t206;
t207 = pkin(7) + t99;
t81 = t207 * t123;
t82 = t207 * t126;
t43 = -t122 * t81 + t125 * t82;
t159 = qJD(4) * t207;
t67 = t123 * t159;
t68 = t126 * t159;
t120 = cos(pkin(9));
t127 = cos(qJ(2));
t172 = qJD(1) * t127;
t124 = sin(qJ(2));
t173 = qJD(1) * t124;
t97 = t118 * t173;
t72 = t120 * t172 - t97;
t86 = t122 * t126 + t123 * t125;
t220 = -qJD(5) * t43 + t122 * t67 - t125 * t68 + t86 * t72;
t143 = t122 * t123 - t125 * t126;
t42 = -t122 * t82 - t125 * t81;
t219 = qJD(5) * t42 - t122 * t68 - t125 * t67 + t143 * t72;
t84 = t118 * t127 + t120 * t124;
t37 = t143 * t84;
t210 = m(6) * pkin(4);
t218 = -mrSges(5,1) - t210;
t171 = qJD(2) * t123;
t161 = mrSges(5,3) * t171;
t94 = qJD(4) * mrSges(5,1) - t161;
t170 = qJD(2) * t126;
t160 = mrSges(5,3) * t170;
t95 = -qJD(4) * mrSges(5,2) + t160;
t144 = -t123 * t94 + t126 * t95;
t217 = -qJD(2) * mrSges(4,2) + t144;
t65 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t89;
t90 = qJDD(2) * t123 + t126 * t164;
t66 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t90;
t216 = -t123 * t66 + t126 * t65;
t110 = t126 * qJD(3);
t168 = qJD(4) * t123;
t182 = qJDD(2) * pkin(2);
t165 = qJD(1) * qJD(2);
t157 = t124 * t165;
t91 = t127 * qJDD(1) - t157;
t80 = t91 + t182;
t158 = t127 * t165;
t92 = qJDD(1) * t124 + t158;
t47 = t118 * t80 + t120 * t92;
t40 = qJDD(2) * pkin(6) + t47;
t96 = qJD(2) * pkin(2) + t172;
t63 = t118 * t96 + t120 * t173;
t59 = qJD(2) * pkin(6) + t63;
t14 = qJD(4) * t110 + t123 * qJDD(3) + t126 * t40 - t168 * t59;
t169 = qJD(3) * t123;
t49 = t126 * t59 + t169;
t15 = -qJD(4) * t49 + t126 * qJDD(3) - t123 * t40;
t215 = -t123 * t15 + t126 * t14;
t115 = qJD(4) + qJD(5);
t76 = t143 * qJD(2);
t77 = t86 * qJD(2);
t45 = mrSges(6,1) * t76 + mrSges(6,2) * t77;
t213 = t225 * qJD(2) + t45;
t211 = qJD(2) ^ 2;
t208 = t77 / 0.2e1;
t205 = pkin(2) * t120;
t155 = pkin(7) * qJD(2) + t59;
t41 = t126 * t155 + t169;
t192 = t122 * t41;
t39 = -t123 * t155 + t110;
t32 = qJD(4) * pkin(4) + t39;
t10 = t125 * t32 - t192;
t200 = t10 * mrSges(6,3);
t199 = t77 * mrSges(6,3);
t198 = t77 * Ifges(6,4);
t106 = cos(t116);
t121 = cos(pkin(8));
t178 = t112 * t121;
t119 = sin(pkin(8));
t179 = t112 * t119;
t180 = t111 * t121;
t181 = t111 * t119;
t197 = (-t106 * t181 - t178) * mrSges(6,1) + (-t106 * t179 + t180) * mrSges(6,2);
t196 = (-t106 * t180 + t179) * mrSges(6,1) + (-t106 * t178 - t181) * mrSges(6,2);
t195 = Ifges(5,4) * t123;
t194 = Ifges(5,4) * t126;
t188 = t125 * t41;
t177 = t119 * t123;
t176 = t119 * t126;
t175 = t121 * t123;
t174 = t121 * t126;
t167 = qJD(4) * t126;
t166 = qJDD(2) * mrSges(4,2);
t163 = pkin(4) * t171;
t162 = pkin(4) * t168;
t46 = -t118 * t92 + t120 * t80;
t62 = t120 * t96 - t97;
t153 = -g(1) * t119 + g(2) * t121;
t152 = mrSges(3,1) * t127 - mrSges(3,2) * t124;
t151 = mrSges(3,1) * t124 + mrSges(3,2) * t127;
t149 = mrSges(5,1) * t123 + mrSges(5,2) * t126;
t148 = -mrSges(6,1) * t111 - mrSges(6,2) * t112;
t147 = t126 * Ifges(5,2) + t195;
t146 = Ifges(5,5) * t126 - Ifges(5,6) * t123;
t11 = t122 * t32 + t188;
t48 = -t123 * t59 + t110;
t145 = -t123 * t48 + t126 * t49;
t83 = t118 * t124 - t120 * t127;
t38 = -qJDD(2) * pkin(3) - t46;
t58 = -qJD(2) * pkin(3) - t62;
t140 = t58 * t149;
t139 = t123 * (Ifges(5,1) * t126 - t195);
t138 = t86 * qJD(5);
t137 = t143 * qJD(5);
t133 = (-t49 * t123 - t48 * t126) * qJD(4) + t215;
t114 = qJDD(4) + qJDD(5);
t8 = qJDD(4) * pkin(4) - pkin(7) * t90 + t15;
t9 = pkin(7) * t89 + t14;
t2 = qJD(5) * t10 + t122 * t8 + t125 * t9;
t27 = -qJD(2) * t137 + t122 * t89 + t125 * t90;
t28 = -qJD(2) * t138 - t122 * t90 + t125 * t89;
t3 = -qJD(5) * t11 - t122 * t9 + t125 * t8;
t34 = -t76 * Ifges(6,2) + t115 * Ifges(6,6) + t198;
t69 = Ifges(6,4) * t76;
t35 = t77 * Ifges(6,1) + t115 * Ifges(6,5) - t69;
t52 = -qJD(2) * t102 - t62;
t132 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t76 * t200 + t34 * t208 - t52 * (mrSges(6,1) * t77 - mrSges(6,2) * t76) + Ifges(6,3) * t114 - t77 * (-Ifges(6,1) * t76 - t198) / 0.2e1 + Ifges(6,6) * t28 + Ifges(6,5) * t27 - t115 * (-Ifges(6,5) * t76 - Ifges(6,6) * t77) / 0.2e1 + (-Ifges(6,2) * t77 + t35 - t69) * t76 / 0.2e1;
t51 = -qJD(4) * t86 - t138;
t103 = Ifges(5,4) * t170;
t93 = -t102 - t205;
t75 = Ifges(5,1) * t171 + Ifges(5,5) * qJD(4) + t103;
t74 = Ifges(5,6) * qJD(4) + qJD(2) * t147;
t73 = t83 * qJD(2);
t71 = t84 * qJD(2);
t61 = mrSges(6,1) * t115 - t199;
t60 = -mrSges(6,2) * t115 - mrSges(6,3) * t76;
t53 = -mrSges(5,1) * t89 + mrSges(5,2) * t90;
t50 = -qJD(4) * t143 - t137;
t36 = t86 * t84;
t24 = -pkin(4) * t89 + t38;
t21 = -mrSges(6,2) * t114 + mrSges(6,3) * t28;
t20 = mrSges(6,1) * t114 - mrSges(6,3) * t27;
t13 = t125 * t39 - t192;
t12 = -t122 * t39 - t188;
t7 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t5 = t115 * t37 + t86 * t73;
t4 = t143 * t73 + t51 * t84;
t1 = [m(2) * qJDD(1) - t36 * t20 - t37 * t21 + t4 * t60 + t5 * t61 + (t53 + t7) * t83 - t151 * t211 - t217 * t73 + t213 * t71 + (-mrSges(4,1) * t83 + t152) * qJDD(2) + (-t166 + (-t123 * t95 - t126 * t94) * qJD(4) + t216) * t84 + (-m(2) - m(3) - t224) * g(3) + m(5) * (t133 * t84 - t145 * t73 + t38 * t83 + t58 * t71) + m(3) * (t124 * t92 + t127 * t91) + m(6) * (t10 * t5 + t11 * t4 - t2 * t37 + t24 * t83 - t3 * t36 + t52 * t71) + m(4) * (-t46 * t83 + t47 * t84 - t62 * t71 - t63 * t73); (-t92 + t158) * mrSges(3,2) + (Ifges(5,1) * t90 + Ifges(5,4) * t221) * t123 + (mrSges(6,2) * t24 - mrSges(6,3) * t3 + Ifges(6,1) * t27 + Ifges(6,4) * t28 + Ifges(6,5) * t114) * t86 + t147 * t221 + (t140 + t146 * qJD(4) / 0.2e1) * qJD(4) + m(4) * (t118 * t47 + t120 * t46) * pkin(2) + qJDD(4) * (Ifges(5,5) * t123 + Ifges(5,6) * t126) + t11 * t51 * mrSges(6,3) + t115 * (Ifges(6,5) * t50 + Ifges(6,6) * t51) / 0.2e1 + t93 * t7 - t76 * (Ifges(6,4) * t50 + Ifges(6,2) * t51) / 0.2e1 + t51 * t34 / 0.2e1 + t52 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + t42 * t20 + t43 * t21 - t47 * mrSges(4,2) + t50 * t35 / 0.2e1 - t74 * t168 / 0.2e1 + t75 * t167 / 0.2e1 - t166 * t206 - t50 * t200 + t126 * (Ifges(5,4) * t90 + Ifges(5,2) * t89) / 0.2e1 + t90 * t194 / 0.2e1 + (t222 * t105 - t223 * t106 - t127 * t226 - t152) * g(3) + (g(1) * t121 + g(2) * t119) * (t151 + t124 * t226 + (-mrSges(5,3) + t222) * t106 + t223 * t105) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (m(5) * t38 + t53) * (-pkin(3) - t205) + (t126 * (-Ifges(5,2) * t123 + t194) + t139) * t164 / 0.2e1 + (-t167 * t48 - t168 * t49 - t201 + t215) * mrSges(5,3) + t219 * t60 + t220 * t61 + (t10 * t220 + t11 * t219 + t52 * t162 + t2 * t43 + t24 * t93 + t3 * t42) * m(6) + (m(5) * t133 - t167 * t94 - t168 * t95 + t216) * t99 + (-m(4) * t63 - m(5) * t145 - t217) * t72 - (-mrSges(6,1) * t24 + mrSges(6,3) * t2 + Ifges(6,4) * t27 + Ifges(6,2) * t28 + Ifges(6,6) * t114) * t143 + t38 * t150 + t45 * t162 + (t91 + t157) * mrSges(3,1) + (Ifges(6,1) * t50 + Ifges(6,4) * t51) * t208 + (t120 * t182 + t46) * mrSges(4,1) + (m(4) * t62 - m(5) * t58 - m(6) * t52 - t213) * t84 * qJD(1); t123 * t65 + t126 * t66 - t143 * t20 + t86 * t21 + t50 * t60 + t51 * t61 + t144 * qJD(4) + (t10 * t51 + t11 * t50 - t143 * t3 + t2 * t86 + t153) * m(6) + (qJD(4) * t145 + t123 * t14 + t126 * t15 + t153) * m(5) + (qJDD(3) + t153) * m(4); t132 + Ifges(5,6) * t89 + Ifges(5,5) * t90 - t13 * t60 - t12 * t61 - t14 * mrSges(5,2) + t15 * mrSges(5,1) + t74 * t171 / 0.2e1 - t146 * t164 / 0.2e1 - t211 * t139 / 0.2e1 - t45 * t163 - m(6) * (t10 * t12 + t11 * t13 + t163 * t52) - qJD(2) * t140 + t11 * t199 + (t122 * t2 + t125 * t3 + (-t10 * t122 + t11 * t125) * qJD(5)) * t210 + Ifges(5,3) * qJDD(4) + (t94 + t161) * t49 + (-t95 + t160) * t48 + (t123 * t210 - t148 + t149) * t201 - (-Ifges(5,2) * t171 + t103 + t75) * t170 / 0.2e1 + (-(-t106 * t176 + t175) * mrSges(5,2) - t197 + t218 * (-t106 * t177 - t174)) * g(2) + (-(-t106 * t174 - t177) * mrSges(5,2) - t196 + t218 * (-t106 * t175 + t176)) * g(1) + ((-t122 * t61 + t125 * t60) * qJD(5) + t122 * t21 + t125 * t20) * pkin(4); t132 + (t61 + t199) * t11 - t148 * t201 - g(2) * t197 - g(1) * t196 - t10 * t60;];
tau = t1;

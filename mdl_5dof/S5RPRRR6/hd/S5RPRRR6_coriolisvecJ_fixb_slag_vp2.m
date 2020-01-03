% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:57
% EndTime: 2019-12-31 19:01:04
% DurationCPUTime: 2.80s
% Computational Cost: add. (3837->324), mult. (9144->455), div. (0->0), fcn. (5885->8), ass. (0->166)
t127 = qJD(3) + qJD(4);
t230 = t127 * Ifges(5,6) / 0.2e1;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t135 = cos(qJ(3));
t120 = sin(pkin(9)) * pkin(1) + pkin(6);
t116 = t120 * qJD(1);
t161 = pkin(7) * qJD(1) + t116;
t132 = sin(qJ(3));
t170 = t132 * qJD(2);
t94 = t135 * t161 + t170;
t138 = qJD(3) * t94;
t182 = t131 * t94;
t126 = t135 * qJD(2);
t146 = t161 * t132;
t93 = t126 - t146;
t90 = qJD(3) * pkin(3) + t93;
t50 = t134 * t90 - t182;
t121 = qJD(3) * t126;
t87 = -qJD(3) * t146 + t121;
t13 = qJD(4) * t50 - t131 * t138 + t134 * t87;
t180 = t134 * t94;
t51 = t131 * t90 + t180;
t49 = pkin(8) * t127 + t51;
t113 = t131 * t132 - t134 * t135;
t108 = t113 * qJD(1);
t114 = t131 * t135 + t132 * t134;
t109 = t114 * qJD(1);
t166 = -cos(pkin(9)) * pkin(1) - pkin(2);
t115 = -pkin(3) * t135 + t166;
t110 = qJD(1) * t115;
t59 = t108 * pkin(4) - t109 * pkin(8) + t110;
t15 = -t130 * t49 + t133 * t59;
t173 = qJD(3) * t132;
t169 = pkin(3) * t173;
t85 = t127 * t113;
t78 = t85 * qJD(1);
t86 = t127 * t114;
t79 = t86 * qJD(1);
t28 = pkin(4) * t79 + pkin(8) * t78 + qJD(1) * t169;
t2 = qJD(5) * t15 + t13 * t133 + t130 * t28;
t16 = t130 * t59 + t133 * t49;
t3 = -qJD(5) * t16 - t13 * t130 + t133 * t28;
t229 = -t130 * t3 + t133 * t2;
t185 = t127 * Ifges(5,5);
t228 = t110 * mrSges(5,2) + t185 / 0.2e1;
t188 = t109 * Ifges(5,4);
t227 = t230 + t188 / 0.2e1 - t108 * Ifges(5,2) / 0.2e1;
t91 = -t109 * t130 + t127 * t133;
t46 = qJD(5) * t91 - t133 * t78;
t92 = t109 * t133 + t127 * t130;
t47 = -qJD(5) * t92 + t130 * t78;
t10 = -mrSges(6,1) * t47 + mrSges(6,2) * t46;
t14 = qJD(4) * t51 + t131 * t87 + t134 * t138;
t224 = m(6) * t14 + t10;
t149 = -t130 * t16 - t133 * t15;
t223 = qJD(5) * t149 + t229;
t105 = qJD(5) + t108;
t156 = mrSges(6,1) * t130 + mrSges(6,2) * t133;
t48 = -pkin(4) * t127 - t50;
t143 = t48 * t156;
t151 = Ifges(6,5) * t133 - Ifges(6,6) * t130;
t193 = Ifges(6,4) * t133;
t153 = -Ifges(6,2) * t130 + t193;
t194 = Ifges(6,4) * t130;
t155 = Ifges(6,1) * t133 - t194;
t206 = t133 / 0.2e1;
t207 = -t130 / 0.2e1;
t212 = t92 / 0.2e1;
t204 = Ifges(6,4) * t92;
t41 = t91 * Ifges(6,2) + t105 * Ifges(6,6) + t204;
t88 = Ifges(6,4) * t91;
t42 = Ifges(6,1) * t92 + Ifges(6,5) * t105 + t88;
t222 = t42 * t206 + t41 * t207 + t105 * t151 / 0.2e1 + t155 * t212 + t91 * t153 / 0.2e1 + t143;
t221 = -t110 * mrSges(5,1) - t15 * mrSges(6,1) + t16 * mrSges(6,2) + t227;
t220 = -t41 / 0.2e1;
t219 = -t130 * t15 + t133 * t16;
t218 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t46 + Ifges(6,6) * t47;
t217 = t46 / 0.2e1;
t216 = t47 / 0.2e1;
t215 = t79 / 0.2e1;
t214 = -t91 / 0.2e1;
t213 = -t92 / 0.2e1;
t210 = -t105 / 0.2e1;
t209 = t108 / 0.2e1;
t208 = -t109 / 0.2e1;
t205 = m(5) * t110;
t198 = pkin(7) + t120;
t111 = t198 * t132;
t112 = t198 * t135;
t145 = -t111 * t134 - t112 * t131;
t201 = t14 * t145;
t200 = t91 * Ifges(6,6);
t199 = t92 * Ifges(6,5);
t190 = t109 * mrSges(5,3);
t197 = mrSges(5,1) * t127 + mrSges(6,1) * t91 - mrSges(6,2) * t92 - t190;
t196 = mrSges(5,3) * t108;
t195 = Ifges(4,4) * t132;
t104 = Ifges(5,4) * t108;
t192 = t105 * Ifges(6,3);
t189 = t109 * Ifges(5,1);
t187 = t113 * t14;
t118 = t166 * qJD(1);
t186 = t118 * mrSges(4,2);
t179 = Ifges(4,5) * qJD(3);
t178 = Ifges(4,6) * qJD(3);
t177 = t108 * t130;
t176 = t108 * t133;
t175 = qJD(1) * t132;
t174 = qJD(1) * t135;
t172 = qJD(5) * t130;
t171 = qJD(5) * t133;
t125 = Ifges(4,4) * t174;
t165 = t179 / 0.2e1;
t164 = -t178 / 0.2e1;
t163 = m(4) * t120 + mrSges(4,3);
t162 = qJD(3) * t198;
t160 = t178 / 0.2e1 + (t135 * Ifges(4,2) + t195) * qJD(1) / 0.2e1 - t118 * mrSges(4,1);
t158 = t135 * t162;
t82 = pkin(4) * t109 + pkin(8) * t108;
t157 = mrSges(6,1) * t133 - mrSges(6,2) * t130;
t154 = Ifges(6,1) * t130 + t193;
t152 = Ifges(6,2) * t133 + t194;
t150 = Ifges(6,5) * t130 + Ifges(6,6) * t133;
t17 = mrSges(6,1) * t79 - mrSges(6,3) * t46;
t18 = -mrSges(6,2) * t79 + mrSges(6,3) * t47;
t147 = -t130 * t17 + t133 * t18;
t70 = pkin(4) * t113 - pkin(8) * t114 + t115;
t77 = -t111 * t131 + t112 * t134;
t29 = -t130 * t77 + t133 * t70;
t30 = t130 * t70 + t133 * t77;
t101 = t116 * t135 + t170;
t60 = -mrSges(6,2) * t105 + mrSges(6,3) * t91;
t61 = mrSges(6,1) * t105 - mrSges(6,3) * t92;
t95 = -mrSges(5,2) * t127 - t196;
t144 = -t130 * t61 + t133 * t60 + t95;
t137 = m(6) * (-t15 * t171 - t16 * t172 + t229) - t61 * t171 - t60 * t172 + t147;
t40 = t192 + t199 + t200;
t68 = -t104 + t185 + t189;
t8 = Ifges(6,4) * t46 + Ifges(6,2) * t47 + Ifges(6,6) * t79;
t9 = t46 * Ifges(6,1) + t47 * Ifges(6,4) + t79 * Ifges(6,5);
t136 = t130 * t9 / 0.2e1 - Ifges(5,5) * t78 - Ifges(5,6) * t79 - t13 * mrSges(5,2) + t177 * t220 - t50 * t196 + t8 * t206 + t150 * t215 + t152 * t216 + t154 * t217 + t42 * t176 / 0.2e1 + (-t104 + t68) * t209 + (-t188 + t40) * t208 + (-mrSges(5,1) - t157) * t14 + (-Ifges(5,1) * t208 - t151 * t210 - t153 * t214 - t155 * t213 + t143 + t228) * t108 + (Ifges(6,5) * t213 - Ifges(5,2) * t209 + Ifges(6,6) * t214 + Ifges(6,3) * t210 + t221 + t230) * t109 + (-t15 * t176 - t16 * t177 + t223) * mrSges(6,3) + t222 * qJD(5);
t119 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t174;
t117 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t175;
t107 = Ifges(4,1) * t175 + t125 + t179;
t103 = t132 * t162;
t100 = -t116 * t132 + t126;
t98 = t101 * qJD(3);
t97 = -t116 * t173 + t121;
t81 = mrSges(5,1) * t108 + mrSges(5,2) * t109;
t73 = Ifges(6,3) * t79;
t66 = pkin(3) * t175 + t82;
t56 = t134 * t93 - t182;
t55 = t131 * t93 + t180;
t45 = pkin(4) * t86 + pkin(8) * t85 + t169;
t34 = qJD(4) * t77 - t103 * t131 + t134 * t158;
t33 = qJD(4) * t145 - t134 * t103 - t131 * t158;
t22 = t130 * t82 + t133 * t50;
t21 = -t130 * t50 + t133 * t82;
t20 = t130 * t66 + t133 * t56;
t19 = -t130 * t56 + t133 * t66;
t5 = -qJD(5) * t30 - t130 * t33 + t133 * t45;
t4 = qJD(5) * t29 + t130 * t45 + t133 * t33;
t1 = [t115 * (mrSges(5,1) * t79 - mrSges(5,2) * t78) + t33 * t95 - t145 * t10 + t4 * t60 + t5 * t61 + t29 * t17 + t30 * t18 + (t40 / 0.2e1 + t192 / 0.2e1 + t200 / 0.2e1 + t199 / 0.2e1 - t221 - t227) * t86 - (-t104 / 0.2e1 + t189 / 0.2e1 + t68 / 0.2e1 + t149 * mrSges(6,3) + t222 + t228) * t85 - t197 * t34 + m(6) * (t15 * t5 + t16 * t4 + t2 * t30 + t29 * t3 + t34 * t48 - t201) + m(5) * (t13 * t77 + t33 * t51 - t34 * t50 - t201) + (t145 * t78 + t50 * t85 - t51 * t86 - t77 * t79) * mrSges(5,3) + (t163 * t97 + (t107 / 0.2e1 - t120 * t117 + 0.3e1 / 0.2e1 * t125 + t165 - t163 * t100 + 0.2e1 * t186) * qJD(3)) * t135 + (t73 / 0.2e1 + Ifges(5,4) * t78 - t13 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t79 + t218) * t113 + (t151 * t215 + t155 * t217 + t153 * t216 - Ifges(5,1) * t78 - Ifges(5,4) * t79 + t9 * t206 + t8 * t207 + (mrSges(5,3) + t156) * t14 + (-t130 * t2 - t133 * t3) * mrSges(6,3) + (-mrSges(6,3) * t219 + t133 * t220 + t150 * t210 + t152 * t214 + t154 * t213 + t48 * t157 + t42 * t207) * qJD(5)) * t114 + (t163 * t98 + (-t120 * t119 + t164 - t163 * t101 + (t166 * mrSges(4,1) - 0.3e1 / 0.2e1 * t195 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t135) * qJD(1) + (0.2e1 * t205 + t81 + qJD(1) * (mrSges(5,1) * t113 + mrSges(5,2) * t114)) * pkin(3) - t160) * qJD(3)) * t132; -t197 * t86 + (-t78 * mrSges(5,3) + t10) * t113 - t144 * t85 + (-t132 * t117 + t135 * t119 + (-t132 ^ 2 - t135 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + (-t79 * mrSges(5,3) + (-t130 * t60 - t133 * t61) * qJD(5) + t147) * t114 + m(4) * (t97 * t132 - t135 * t98 + (-t100 * t132 + t101 * t135) * qJD(3)) + m(5) * (t114 * t13 - t50 * t86 - t51 * t85 + t187) + m(6) * (t114 * t223 - t219 * t85 + t48 * t86 + t187); t101 * t117 - t100 * t119 - t97 * mrSges(4,2) - t98 * mrSges(4,1) - t56 * t95 + t51 * t190 - t20 * t60 - t19 * t61 + t136 + t197 * t55 - m(6) * (t15 * t19 + t16 * t20 + t48 * t55) - m(5) * (-t50 * t55 + t51 * t56) + ((-t186 + t165 - t125 / 0.2e1 - t107 / 0.2e1 + t100 * mrSges(4,3)) * t135 + (t164 + t101 * mrSges(4,3) + (t195 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t135) * qJD(1) + (-t81 - t205) * pkin(3) + t160) * t132) * qJD(1) + (m(5) * (t13 * t131 - t134 * t14) + (-t131 * t79 + t134 * t78) * mrSges(5,3) + ((-m(5) * t50 + m(6) * t48 - t197) * t131 + (m(5) * t51 + m(6) * t219 + t144) * t134) * qJD(4)) * pkin(3) + t137 * (pkin(3) * t131 + pkin(8)) + t224 * (-pkin(3) * t134 - pkin(4)); -t50 * t95 - t22 * t60 - t21 * t61 + t136 + (t190 + t197) * t51 - m(6) * (t15 * t21 + t16 * t22 + t48 * t51) + t137 * pkin(8) - t224 * pkin(4); t73 - t48 * (mrSges(6,1) * t92 + mrSges(6,2) * t91) + (Ifges(6,1) * t91 - t204) * t213 + t41 * t212 + (Ifges(6,5) * t91 - Ifges(6,6) * t92) * t210 - t15 * t60 + t16 * t61 + (t15 * t91 + t16 * t92) * mrSges(6,3) + (-Ifges(6,2) * t92 + t42 + t88) * t214 + t218;];
tauc = t1(:);

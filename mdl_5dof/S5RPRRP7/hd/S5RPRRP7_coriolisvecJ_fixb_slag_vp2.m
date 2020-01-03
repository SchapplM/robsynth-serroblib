% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:34
% EndTime: 2019-12-31 18:44:41
% DurationCPUTime: 3.51s
% Computational Cost: add. (1996->349), mult. (4905->448), div. (0->0), fcn. (2656->6), ass. (0->166)
t211 = Ifges(5,1) + Ifges(6,1);
t206 = Ifges(5,5) + Ifges(6,4);
t212 = qJD(3) / 0.2e1;
t210 = Ifges(5,6) - Ifges(6,6);
t143 = qJD(1) * qJD(3);
t98 = sin(qJ(3));
t133 = t98 * t143;
t142 = qJD(3) * qJD(4);
t100 = cos(qJ(3));
t144 = qJD(3) * t100;
t97 = sin(qJ(4));
t150 = qJD(4) * t97;
t99 = cos(qJ(4));
t51 = t99 * t142 + (t144 * t99 - t150 * t98) * qJD(1);
t149 = qJD(4) * t99;
t52 = t97 * t142 + (t144 * t97 + t149 * t98) * qJD(1);
t209 = (-Ifges(5,4) + Ifges(6,5)) * t52 + t211 * t51 + t206 * t133;
t170 = Ifges(6,5) * t97;
t175 = Ifges(5,4) * t97;
t208 = t211 * t99 + t170 - t175;
t153 = qJD(1) * t98;
t207 = t153 / 0.2e1;
t132 = Ifges(4,5) * t212;
t92 = sin(pkin(8)) * pkin(1) + pkin(6);
t85 = t92 * qJD(1);
t65 = qJD(2) * t100 - t98 * t85;
t61 = qJD(3) * t65;
t148 = t98 * qJD(2);
t66 = t100 * t85 + t148;
t60 = qJD(3) * pkin(7) + t66;
t137 = -cos(pkin(8)) * pkin(1) - pkin(2);
t75 = -pkin(3) * t100 - t98 * pkin(7) + t137;
t63 = t75 * qJD(1);
t20 = -t60 * t97 + t63 * t99;
t21 = t60 * t99 + t63 * t97;
t110 = t20 * t99 + t21 * t97;
t200 = qJD(5) - t20;
t146 = qJD(1) * t100;
t91 = qJD(4) - t146;
t14 = -pkin(4) * t91 + t200;
t15 = qJ(5) * t91 + t21;
t111 = t14 * t99 - t15 * t97;
t169 = Ifges(6,5) * t99;
t115 = Ifges(6,3) * t97 + t169;
t174 = Ifges(5,4) * t99;
t119 = -Ifges(5,2) * t97 + t174;
t124 = mrSges(6,1) * t97 - mrSges(6,3) * t99;
t126 = mrSges(5,1) * t97 + mrSges(5,2) * t99;
t167 = Ifges(6,6) * t97;
t168 = Ifges(5,6) * t97;
t172 = Ifges(5,5) * t99;
t173 = Ifges(6,4) * t99;
t183 = t99 / 0.2e1;
t185 = t97 / 0.2e1;
t186 = -t97 / 0.2e1;
t152 = qJD(3) * t97;
t80 = t153 * t99 + t152;
t189 = t80 / 0.2e1;
t59 = -qJD(3) * pkin(3) - t65;
t147 = t99 * qJD(3);
t79 = t153 * t97 - t147;
t19 = t79 * pkin(4) - t80 * qJ(5) + t59;
t191 = t79 / 0.2e1;
t192 = -t79 / 0.2e1;
t171 = Ifges(6,5) * t79;
t77 = Ifges(5,4) * t79;
t202 = t206 * t91 + t211 * t80 + t171 - t77;
t203 = t91 / 0.2e1;
t76 = Ifges(6,5) * t80;
t24 = Ifges(6,6) * t91 + Ifges(6,3) * t79 + t76;
t176 = Ifges(5,4) * t80;
t27 = -Ifges(5,2) * t79 + Ifges(5,6) * t91 + t176;
t101 = t111 * mrSges(6,2) - t110 * mrSges(5,3) + t115 * t191 + t119 * t192 + t19 * t124 + t59 * t126 + t24 * t185 + t27 * t186 + t208 * t189 + (t167 + t173 - t168 + t172) * t203 + t202 * t183;
t93 = Ifges(4,4) * t146;
t205 = t101 - t65 * mrSges(4,3) + Ifges(4,1) * t207 + t132 + t93 / 0.2e1;
t204 = t211 * t97 - t169 + t174;
t131 = -Ifges(4,6) * qJD(3) / 0.2e1;
t156 = t100 * t92;
t201 = t99 * t156 + t97 * t75;
t199 = t206 * t97 + t210 * t99;
t128 = pkin(3) * t98 - pkin(7) * t100;
t83 = t128 * qJD(3);
t74 = qJD(1) * t83;
t3 = t63 * t149 - t150 * t60 + t99 * t61 + t97 * t74;
t4 = -qJD(4) * t21 - t61 * t97 + t74 * t99;
t129 = t3 * t99 - t4 * t97;
t1 = qJ(5) * t133 + qJD(5) * t91 + t3;
t2 = -pkin(4) * t133 - t4;
t130 = t1 * t99 + t2 * t97;
t197 = qJD(4) * t201 - t99 * t83;
t138 = -Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t139 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t140 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t87 = t137 * qJD(1);
t196 = t138 * t79 - t139 * t91 - t140 * t80 - t15 * mrSges(6,3) - t20 * mrSges(5,1) - t87 * mrSges(4,1) - Ifges(5,6) * t192 - Ifges(6,6) * t191 - t131 + (Ifges(4,4) * t98 + t100 * Ifges(4,2)) * qJD(1) / 0.2e1 + t14 * mrSges(6,1) + t21 * mrSges(5,2) - (Ifges(5,3) + Ifges(6,2)) * t203 - t206 * t189;
t195 = t51 / 0.2e1;
t194 = -t52 / 0.2e1;
t193 = t52 / 0.2e1;
t190 = -t80 / 0.2e1;
t188 = -t91 / 0.2e1;
t184 = -t99 / 0.2e1;
t178 = mrSges(5,3) * t79;
t177 = mrSges(5,3) * t80;
t165 = t87 * mrSges(4,2);
t164 = t92 * t97;
t163 = t99 * t75;
t54 = -mrSges(5,2) * t91 - t178;
t57 = -t79 * mrSges(6,2) + mrSges(6,3) * t91;
t160 = t54 + t57;
t55 = mrSges(5,1) * t91 - t177;
t56 = -mrSges(6,1) * t91 + mrSges(6,2) * t80;
t159 = -t55 + t56;
t82 = t128 * qJD(1);
t31 = t99 * t65 + t97 * t82;
t158 = t75 * t149 + t97 * t83;
t141 = mrSges(4,3) * t153;
t157 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t79 - mrSges(5,2) * t80 - t141;
t62 = qJD(3) * t66;
t151 = qJD(3) * t98;
t136 = mrSges(4,3) * t146;
t135 = m(4) * t92 + mrSges(4,3);
t134 = pkin(4) + t164;
t127 = mrSges(5,1) * t99 - mrSges(5,2) * t97;
t125 = mrSges(6,1) * t99 + mrSges(6,3) * t97;
t118 = Ifges(5,2) * t99 + t175;
t114 = -Ifges(6,3) * t99 + t170;
t113 = pkin(4) * t99 + qJ(5) * t97;
t112 = pkin(4) * t97 - qJ(5) * t99;
t30 = -t65 * t97 + t82 * t99;
t37 = -mrSges(6,1) * t133 + t51 * mrSges(6,2);
t108 = t112 + t92;
t35 = -mrSges(6,2) * t52 + mrSges(6,3) * t133;
t36 = mrSges(5,1) * t133 - mrSges(5,3) * t51;
t38 = -mrSges(5,2) * t133 - mrSges(5,3) * t52;
t106 = (t35 + t38) * t99 + (-t36 + t37) * t97;
t105 = t159 * t99 - t160 * t97;
t104 = -t4 * mrSges(5,1) + t2 * mrSges(6,1) + t3 * mrSges(5,2) - t1 * mrSges(6,3);
t90 = Ifges(6,2) * t133;
t89 = Ifges(5,3) * t133;
t88 = -qJD(3) * mrSges(4,2) + t136;
t84 = -pkin(3) - t113;
t69 = qJD(4) * t112 - qJD(5) * t97;
t53 = t108 * t98;
t48 = Ifges(6,4) * t51;
t47 = Ifges(5,5) * t51;
t46 = Ifges(5,6) * t52;
t45 = Ifges(6,6) * t52;
t42 = mrSges(6,1) * t79 - mrSges(6,3) * t80;
t41 = pkin(4) * t80 + qJ(5) * t79;
t39 = -t156 * t97 + t163;
t34 = t100 * t134 - t163;
t33 = -qJ(5) * t100 + t201;
t32 = t148 + (qJD(1) * t112 + t85) * t100;
t23 = -pkin(4) * t153 - t30;
t22 = qJ(5) * t153 + t31;
t18 = (qJD(4) * t113 - qJD(5) * t99) * t98 + t108 * t144;
t17 = mrSges(5,1) * t52 + mrSges(5,2) * t51;
t16 = mrSges(6,1) * t52 - mrSges(6,3) * t51;
t11 = t51 * Ifges(5,4) - t52 * Ifges(5,2) + Ifges(5,6) * t133;
t10 = t51 * Ifges(6,5) + Ifges(6,6) * t133 + t52 * Ifges(6,3);
t9 = t151 * t164 - t197;
t8 = (-t100 * t150 - t147 * t98) * t92 + t158;
t7 = -t134 * t151 + t197;
t6 = (-t150 * t92 - qJD(5)) * t100 + (-t92 * t99 + qJ(5)) * t151 + t158;
t5 = t52 * pkin(4) - t51 * qJ(5) - t80 * qJD(5) + t62;
t12 = [t53 * t16 + t18 * t42 + t33 * t35 + t34 * t37 + t39 * t36 + t201 * t38 + t8 * t54 + t9 * t55 + t7 * t56 + t6 * t57 + m(5) * (t20 * t9 + t201 * t3 + t21 * t8 + t4 * t39) + m(6) * (t1 * t33 + t14 * t7 + t15 * t6 + t18 * t19 + t2 * t34 + t5 * t53) + (-t89 / 0.2e1 - t90 / 0.2e1 - t48 / 0.2e1 - t45 / 0.2e1 + t46 / 0.2e1 - t47 / 0.2e1 + t135 * t61 + t138 * t52 - t140 * t51 + (0.3e1 / 0.2e1 * t93 + 0.2e1 * t165 + t132 + t205) * qJD(3) + t104) * t100 + (t119 * t194 + t115 * t193 + t10 * t185 + t11 * t186 + t5 * t124 + (-t3 * t97 - t4 * t99) * mrSges(5,3) + (-t1 * t97 + t2 * t99) * mrSges(6,2) + (mrSges(4,3) + t126) * t62 + (-t135 * t66 + t131 - t196) * qJD(3) + (t118 * t191 + t114 * t192 + t59 * t127 + t19 * t125 + t27 * t184 + (t20 * t97 - t21 * t99) * mrSges(5,3) + (-t14 * t97 - t15 * t99) * mrSges(6,2) + t204 * t190 + t199 * t188 + t202 * t186) * qJD(4) + (t137 * mrSges(4,1) + (-0.3e1 / 0.2e1 * Ifges(4,4) + t173 / 0.2e1 + t167 / 0.2e1 + t172 / 0.2e1 - t168 / 0.2e1) * t98 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t139) * t100) * t143 + t208 * t195 + (qJD(4) * t24 + t209) * t183) * t98 + ((-m(4) * t65 + m(5) * t59 - t157) * t144 + (t17 + (m(5) + m(4)) * t62 - t88 * qJD(3)) * t98) * t92; (-t16 - t17 + (t159 * t97 + t160 * t99 - t136 + t88) * qJD(3) + m(5) * (t147 * t21 - t152 * t20 - t62) + m(6) * (t14 * t152 + t147 * t15 - t5)) * t100 + (t105 * qJD(4) + (t42 - t141 - t157) * qJD(3) + m(5) * (qJD(3) * t59 - t149 * t20 - t150 * t21 + t129) + m(6) * (qJD(3) * t19 + t14 * t149 - t15 * t150 + t130) + t106) * t98; t204 * t195 - m(5) * (t20 * t30 + t21 * t31 + t59 * t66) - m(6) * (t14 * t23 + t15 * t22 + t19 * t32) + (t69 - t32) * t42 + m(6) * (t130 * pkin(7) + t19 * t69 + t5 * t84) + m(5) * (-pkin(3) * t62 + t129 * pkin(7)) + t129 * mrSges(5,3) + t130 * mrSges(6,2) + (-mrSges(4,1) - t127) * t62 - t5 * t125 + ((-m(5) * t110 + m(6) * t111 + t105) * pkin(7) + t101) * qJD(4) + t84 * t16 - t65 * t88 + t157 * t66 - t31 * t54 - t30 * t55 - t23 * t56 - t22 * t57 - t61 * mrSges(4,2) - pkin(3) * t17 + ((t66 * mrSges(4,3) + Ifges(4,4) * t207 + t199 * t212 + t131 + t196) * t98 + (-t165 - t93 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t153 + t132 - t205) * t100) * qJD(1) + t106 * pkin(7) + t209 * t185 + t114 * t193 + t118 * t194 + t11 * t183 + t10 * t184; (-t159 + t177) * t21 + (t14 * t79 + t15 * t80) * mrSges(6,2) - t19 * (t80 * mrSges(6,1) + t79 * mrSges(6,3)) - t59 * (mrSges(5,1) * t80 - mrSges(5,2) * t79) + t89 + t90 + t48 + t45 - t46 + t47 - t104 + qJD(5) * t57 + qJ(5) * t35 - pkin(4) * t37 - t41 * t42 + (-t160 - t178) * t20 + t27 * t189 + (Ifges(6,3) * t80 - t171) * t192 + (-t206 * t79 - t210 * t80) * t188 + (-pkin(4) * t2 + qJ(5) * t1 - t14 * t21 + t200 * t15 - t19 * t41) * m(6) + (-Ifges(5,2) * t80 + t202 - t77) * t191 + (-t211 * t79 - t176 + t24 + t76) * t190; t80 * t42 - t91 * t57 + 0.2e1 * (t2 / 0.2e1 + t15 * t188 + t19 * t189) * m(6) + t37;];
tauc = t12(:);

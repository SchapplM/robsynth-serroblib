% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:40
% DurationCPUTime: 3.26s
% Computational Cost: add. (1938->342), mult. (4355->443), div. (0->0), fcn. (2210->4), ass. (0->169)
t219 = Ifges(5,1) + Ifges(6,1);
t215 = Ifges(6,4) + Ifges(5,5);
t221 = Ifges(5,6) - Ifges(6,6);
t92 = sin(qJ(3));
t94 = cos(qJ(3));
t79 = pkin(3) * t92 - pkin(7) * t94 + qJ(2);
t64 = t79 * qJD(1);
t95 = -pkin(1) - pkin(6);
t86 = qJD(1) * t95 + qJD(2);
t78 = t92 * t86;
t66 = qJD(3) * pkin(7) + t78;
t91 = sin(qJ(4));
t93 = cos(qJ(4));
t26 = t64 * t93 - t66 * t91;
t27 = t64 * t91 + t66 * t93;
t109 = t26 * t93 + t27 * t91;
t208 = qJD(5) - t26;
t153 = qJD(1) * t92;
t89 = qJD(4) + t153;
t17 = -pkin(4) * t89 + t208;
t18 = qJ(5) * t89 + t27;
t111 = t17 * t93 - t18 * t91;
t171 = Ifges(6,5) * t93;
t115 = Ifges(6,3) * t91 + t171;
t176 = Ifges(5,4) * t93;
t121 = -Ifges(5,2) * t91 + t176;
t126 = mrSges(6,1) * t91 - mrSges(6,3) * t93;
t128 = mrSges(5,1) * t91 + mrSges(5,2) * t93;
t169 = Ifges(6,6) * t91;
t170 = Ifges(5,6) * t91;
t174 = Ifges(5,5) * t93;
t175 = Ifges(6,4) * t93;
t187 = t93 / 0.2e1;
t189 = t91 / 0.2e1;
t163 = t86 * t94;
t67 = -qJD(3) * pkin(3) - t163;
t146 = t93 * qJD(3);
t152 = qJD(1) * t94;
t74 = t152 * t91 - t146;
t75 = qJD(3) * t91 + t152 * t93;
t19 = pkin(4) * t74 - qJ(5) * t75 + t67;
t190 = -t91 / 0.2e1;
t191 = t89 / 0.2e1;
t193 = t75 / 0.2e1;
t195 = t74 / 0.2e1;
t196 = -t74 / 0.2e1;
t69 = Ifges(6,5) * t75;
t20 = Ifges(6,6) * t89 + Ifges(6,3) * t74 + t69;
t172 = Ifges(6,5) * t91;
t177 = Ifges(5,4) * t91;
t206 = t219 * t93 + t172 - t177;
t173 = Ifges(6,5) * t74;
t70 = Ifges(5,4) * t74;
t210 = t215 * t89 + t219 * t75 + t173 - t70;
t178 = Ifges(5,4) * t75;
t23 = -Ifges(5,2) * t74 + Ifges(5,6) * t89 + t178;
t220 = t210 * t187 + t189 * t20 + t190 * t23 + t115 * t195 + t121 * t196 + t19 * t126 + t67 * t128 + (t169 + t175 - t170 + t174) * t191 + t206 * t193 + t111 * mrSges(6,2) - t109 * mrSges(5,3);
t155 = Ifges(4,5) * qJD(3);
t180 = Ifges(4,4) * t92;
t216 = qJD(1) / 0.2e1;
t218 = t155 / 0.2e1 + (t94 * Ifges(4,1) - t180) * t216 + t220;
t145 = qJD(1) * qJD(3);
t135 = t94 * t145;
t148 = qJD(4) * t91;
t47 = qJD(4) * t146 + (-t146 * t92 - t148 * t94) * qJD(1);
t48 = -t145 * t91 * t92 + qJD(4) * t75;
t217 = (-Ifges(5,4) + Ifges(6,5)) * t48 + t219 * t47 + t215 * t135;
t112 = pkin(4) * t91 - qJ(5) * t93;
t214 = -qJD(5) * t91 + t89 * t112 - t78;
t213 = t215 * t91 + t221 * t93;
t212 = t219 * t91 - t171 + t176;
t211 = (qJ(2) * (m(3) + m(4)) + mrSges(3,3)) * qJD(1);
t133 = -Ifges(4,6) * qJD(3) / 0.2e1;
t162 = t92 * t95;
t209 = t93 * t162 + t91 * t79;
t149 = qJD(3) * t94;
t142 = t86 * t149;
t147 = qJD(4) * t93;
t132 = pkin(3) * t94 + pkin(7) * t92;
t73 = qJD(3) * t132 + qJD(2);
t58 = t73 * qJD(1);
t4 = t93 * t142 + t64 * t147 - t148 * t66 + t91 * t58;
t204 = qJD(4) * t27 - t58 * t93;
t5 = -t142 * t91 - t204;
t130 = t4 * t93 - t5 * t91;
t1 = qJ(5) * t135 + qJD(5) * t89 + t4;
t103 = (-pkin(4) * qJD(1) + t86 * t91) * t94;
t2 = qJD(3) * t103 + t204;
t131 = t1 * t93 + t2 * t91;
t203 = qJD(4) * t209 - t73 * t93;
t181 = mrSges(5,3) * t75;
t50 = mrSges(5,1) * t89 - t181;
t51 = -mrSges(6,1) * t89 + mrSges(6,2) * t75;
t159 = -t50 + t51;
t166 = t74 * mrSges(5,3);
t49 = -mrSges(5,2) * t89 - t166;
t52 = -mrSges(6,2) * t74 + mrSges(6,3) * t89;
t160 = t49 + t52;
t101 = t159 * t93 - t160 * t91;
t202 = -m(5) * t109 + m(6) * t111 + t101;
t137 = -Ifges(5,3) / 0.2e1 - Ifges(6,2) / 0.2e1;
t138 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t139 = -Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t179 = Ifges(4,4) * t94;
t201 = t137 * t89 - t138 * t74 + t139 * t75 - t18 * mrSges(6,3) - t26 * mrSges(5,1) - Ifges(5,6) * t196 - Ifges(6,6) * t195 - t133 + (-t92 * Ifges(4,2) + t179) * t216 + t17 * mrSges(6,1) + t27 * mrSges(5,2) - t215 * t193 - (Ifges(5,3) + Ifges(6,2)) * t191;
t199 = t47 / 0.2e1;
t198 = -t48 / 0.2e1;
t197 = t48 / 0.2e1;
t194 = -t75 / 0.2e1;
t192 = -t89 / 0.2e1;
t188 = -t93 / 0.2e1;
t77 = t132 * qJD(1);
t165 = t77 * t93;
t164 = t79 * t93;
t35 = t93 * t163 + t91 * t77;
t158 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t74 - mrSges(5,2) * t75 - mrSges(4,3) * t152;
t157 = qJ(2) * mrSges(4,1);
t156 = qJ(2) * mrSges(4,2);
t151 = qJD(3) * mrSges(4,2);
t150 = qJD(3) * t92;
t141 = t95 * t149;
t144 = t93 * t141 + t79 * t147 + t91 * t73;
t143 = t86 * t150;
t140 = t95 * t148;
t136 = t91 * t95 - pkin(4);
t134 = -t155 / 0.2e1;
t129 = mrSges(5,1) * t93 - mrSges(5,2) * t91;
t127 = mrSges(6,1) * t93 + mrSges(6,3) * t91;
t120 = Ifges(5,2) * t93 + t177;
t114 = -Ifges(6,3) * t93 + t172;
t113 = pkin(4) * t93 + qJ(5) * t91;
t110 = t17 * t91 + t18 * t93;
t108 = t26 * t91 - t27 * t93;
t32 = -mrSges(6,1) * t135 + t47 * mrSges(6,2);
t106 = t112 - t95;
t30 = -mrSges(6,2) * t48 + mrSges(6,3) * t135;
t31 = mrSges(5,1) * t135 - mrSges(5,3) * t47;
t33 = -mrSges(5,2) * t135 - mrSges(5,3) * t48;
t102 = (t30 + t33) * t93 + (-t31 + t32) * t91;
t100 = t5 * mrSges(5,1) - t2 * mrSges(6,1) - t4 * mrSges(5,2) + t1 * mrSges(6,3);
t88 = Ifges(6,2) * t135;
t87 = Ifges(5,3) * t135;
t82 = -mrSges(4,3) * t153 - t151;
t80 = -pkin(3) - t113;
t76 = qJD(1) * (t92 * mrSges(4,1) + t94 * mrSges(4,2));
t56 = t106 * t94;
t53 = -t162 * t91 + t164;
t46 = Ifges(6,4) * t47;
t45 = Ifges(5,5) * t47;
t44 = Ifges(5,6) * t48;
t43 = Ifges(6,6) * t48;
t41 = t136 * t92 - t164;
t40 = qJ(5) * t92 + t209;
t37 = mrSges(6,1) * t74 - mrSges(6,3) * t75;
t36 = pkin(4) * t75 + qJ(5) * t74;
t34 = -t163 * t91 + t165;
t29 = t103 - t165;
t28 = qJ(5) * t152 + t35;
t16 = (qJD(4) * t113 - qJD(5) * t93) * t94 - t106 * t150;
t15 = -t141 * t91 - t203;
t14 = -t140 * t92 + t144;
t13 = mrSges(5,1) * t48 + mrSges(5,2) * t47;
t12 = mrSges(6,1) * t48 - mrSges(6,3) * t47;
t11 = t136 * t149 + t203;
t8 = t47 * Ifges(5,4) - t48 * Ifges(5,2) + Ifges(5,6) * t135;
t7 = t47 * Ifges(6,5) + Ifges(6,6) * t135 + t48 * Ifges(6,3);
t6 = qJ(5) * t149 + (qJD(5) - t140) * t92 + t144;
t3 = pkin(4) * t48 - qJ(5) * t47 - qJD(5) * t75 + t143;
t9 = [t11 * t51 + t56 * t12 + t14 * t49 + t15 * t50 + t16 * t37 + t40 * t30 + t53 * t31 + t41 * t32 + t209 * t33 + t6 * t52 + m(5) * (t14 * t27 + t15 * t26 + t209 * t4 + t5 * t53) + m(6) * (t1 * t40 + t11 * t17 + t16 * t19 + t18 * t6 + t2 * t41 + t3 * t56) + (t76 + 0.2e1 * t211) * qJD(2) + (t87 / 0.2e1 + t88 / 0.2e1 + t46 / 0.2e1 - t44 / 0.2e1 + t45 / 0.2e1 + t43 / 0.2e1 + qJD(1) * qJD(2) * mrSges(4,1) + t138 * t48 - t139 * t47 + ((0.3e1 / 0.2e1 * t180 - 0.2e1 * t156) * qJD(1) + (m(5) * t67 - t158) * t95 + t134 - t218) * qJD(3) + t100) * t92 + (t7 * t189 - t95 * t13 + t115 * t197 + t121 * t198 + t3 * t126 + (-t4 * t91 - t5 * t93) * mrSges(5,3) + (-t1 * t91 + t2 * t93) * mrSges(6,2) + (t95 * t82 + t133 + (-m(5) * t95 + t128) * t78 - t201) * qJD(3) + (-t110 * mrSges(6,2) + mrSges(5,3) * t108 + t114 * t196 + t120 * t195 + t127 * t19 + t129 * t67 + t188 * t23 + t213 * t192 + t212 * t194) * qJD(4) + (qJD(2) * mrSges(4,2) + (0.2e1 * t157 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t175 / 0.2e1 + t169 / 0.2e1 + t174 / 0.2e1 - t170 / 0.2e1) * t94 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) - t137) * t92) * qJD(3)) * qJD(1) + t206 * t199 + (t210 * qJD(4) + t8) * t190 + (t20 * qJD(4) + t217) * t187) * t94; (-t12 - t13 - m(6) * t3 + (-m(5) * t108 + m(6) * t110 + t159 * t91 + t160 * t93 + t82) * qJD(3)) * t94 + ((t37 - t158) * qJD(3) + t101 * qJD(4) + m(5) * (qJD(3) * t67 - t147 * t26 - t148 * t27 + t130 - t142) + m(6) * (qJD(3) * t19 + t147 * t17 - t148 * t18 + t131) + t102) * t92 + (t202 - t76 - t211) * qJD(1); ((t133 + (-t157 + t179 / 0.2e1) * qJD(1) + t213 * qJD(3) / 0.2e1 + t201) * t94 + ((t156 - t180 / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t94) * qJD(1) + t134 + t218) * t92) * qJD(1) + (t202 * pkin(7) + t220) * qJD(4) + ((-t82 - t151) * t94 + ((-mrSges(4,1) - t129) * qJD(3) + t158) * t92) * t86 + t214 * t37 + t80 * t12 + t130 * mrSges(5,3) + t131 * mrSges(6,2) - t3 * t127 - t29 * t51 - t28 * t52 - t35 * t49 - t34 * t50 - pkin(3) * t13 + t102 * pkin(7) + t8 * t187 + t7 * t188 + t114 * t197 + t120 * t198 + t212 * t199 + t217 * t189 + (t131 * pkin(7) - t17 * t29 - t18 * t28 + t214 * t19 + t3 * t80) * m(6) + (-pkin(3) * t143 + t130 * pkin(7) - t26 * t34 - t27 * t35 - t67 * t78) * m(5); (t17 * t74 + t18 * t75) * mrSges(6,2) + (-t159 + t181) * t27 + (-t160 - t166) * t26 + t87 + t88 + t46 - t44 + t45 + t43 - t19 * (mrSges(6,1) * t75 + mrSges(6,3) * t74) - t67 * (mrSges(5,1) * t75 - mrSges(5,2) * t74) + qJD(5) * t52 + qJ(5) * t30 - pkin(4) * t32 - t36 * t37 + t100 + t23 * t193 + (Ifges(6,3) * t75 - t173) * t196 + (-t215 * t74 - t221 * t75) * t192 + (-t2 * pkin(4) + t1 * qJ(5) - t17 * t27 + t208 * t18 - t19 * t36) * m(6) + (-Ifges(5,2) * t75 + t210 - t70) * t195 + (-t219 * t74 - t178 + t20 + t69) * t194; t75 * t37 - t89 * t52 + 0.2e1 * (t2 / 0.2e1 + t18 * t192 + t19 * t193) * m(6) + t32;];
tauc = t9(:);

% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:09
% EndTime: 2019-12-31 18:42:19
% DurationCPUTime: 4.19s
% Computational Cost: add. (1995->340), mult. (4949->449), div. (0->0), fcn. (2722->6), ass. (0->157)
t221 = Ifges(5,4) + Ifges(6,4);
t220 = Ifges(5,1) + Ifges(6,1);
t208 = Ifges(5,5) + Ifges(6,5);
t219 = Ifges(5,2) + Ifges(6,2);
t207 = Ifges(6,6) + Ifges(5,6);
t107 = cos(qJ(4));
t218 = t221 * t107;
t105 = sin(qJ(4));
t217 = t221 * t105;
t216 = qJD(3) / 0.2e1;
t106 = sin(qJ(3));
t145 = qJD(1) * qJD(3);
t135 = t106 * t145;
t144 = qJD(3) * qJD(4);
t149 = qJD(4) * t105;
t108 = cos(qJ(3));
t150 = qJD(3) * t108;
t51 = t107 * t144 + (-t106 * t149 + t107 * t150) * qJD(1);
t148 = qJD(4) * t107;
t113 = t105 * t150 + t106 * t148;
t52 = -qJD(1) * t113 - t105 * t144;
t215 = t208 * t135 + t220 * t51 + t221 * t52;
t151 = qJD(3) * t107;
t155 = qJD(1) * t106;
t86 = -t105 * t155 + t151;
t214 = t221 * t86;
t213 = -t219 * t105 + t218;
t212 = t220 * t107 - t217;
t152 = qJD(3) * t105;
t87 = t107 * t155 + t152;
t211 = t221 * t87;
t154 = qJD(1) * t108;
t102 = Ifges(4,4) * t154;
t137 = Ifges(4,5) * t216;
t146 = t106 * qJD(2);
t100 = sin(pkin(8)) * pkin(1) + pkin(6);
t91 = t100 * qJD(1);
t66 = t108 * t91 + t146;
t60 = qJD(3) * pkin(7) + t66;
t138 = -cos(pkin(8)) * pkin(1) - pkin(2);
t81 = -pkin(3) * t108 - t106 * pkin(7) + t138;
t63 = t81 * qJD(1);
t19 = -t105 * t60 + t107 * t63;
t20 = t105 * t63 + t107 * t60;
t115 = t20 * t105 + t19 * t107;
t126 = mrSges(6,1) * t105 + mrSges(6,2) * t107;
t128 = mrSges(5,1) * t105 + mrSges(5,2) * t107;
t162 = Ifges(6,6) * t105;
t163 = Ifges(5,6) * t105;
t164 = Ifges(6,5) * t107;
t165 = Ifges(5,5) * t107;
t184 = t107 / 0.2e1;
t187 = -t105 / 0.2e1;
t190 = t87 / 0.2e1;
t99 = qJD(4) - t154;
t199 = t208 * t99 + t220 * t87 + t214;
t200 = t99 / 0.2e1;
t201 = t86 / 0.2e1;
t205 = t207 * t99 + t219 * t86 + t211;
t65 = qJD(2) * t108 - t106 * t91;
t59 = -qJD(3) * pkin(3) - t65;
t33 = -t86 * pkin(4) + qJD(5) + t59;
t8 = -qJ(5) * t87 + t19;
t7 = pkin(4) * t99 + t8;
t9 = qJ(5) * t86 + t20;
t197 = t115 * mrSges(5,3) + (t9 * t105 + t7 * t107) * mrSges(6,3) - t33 * t126 - t59 * t128 - t213 * t201 - t212 * t190 - (-t162 + t164 - t163 + t165) * t200 - t205 * t187 - t199 * t184;
t209 = t155 / 0.2e1;
t210 = -t65 * mrSges(4,3) + Ifges(4,1) * t209 + t102 / 0.2e1 + t137 - t197;
t206 = t207 * t135 + t219 * t52 + t221 * t51;
t61 = qJD(3) * t65;
t204 = t208 * t105 + t207 * t107;
t203 = t219 * t107 + t217;
t202 = t220 * t105 + t218;
t136 = -Ifges(4,6) * qJD(3) / 0.2e1;
t133 = pkin(3) * t106 - pkin(7) * t108;
t90 = t133 * qJD(3);
t80 = qJD(1) * t90;
t3 = t105 * t80 + t107 * t61 + t63 * t148 - t149 * t60;
t4 = -qJD(4) * t20 - t105 * t61 + t107 * t80;
t131 = -t105 * t4 + t107 * t3;
t141 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t142 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t143 = Ifges(5,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t93 = t138 * qJD(1);
t196 = -t141 * t99 - t142 * t86 - t143 * t87 - t19 * mrSges(5,1) - t7 * mrSges(6,1) - t93 * mrSges(4,1) - t136 + (Ifges(4,4) * t106 + t108 * Ifges(4,2)) * qJD(1) / 0.2e1 + t20 * mrSges(5,2) + t66 * mrSges(4,3) + t9 * mrSges(6,2) - t207 * t201 - (Ifges(6,3) + Ifges(5,3)) * t200 - t208 * t190;
t195 = t51 / 0.2e1;
t194 = t52 / 0.2e1;
t193 = -t86 / 0.2e1;
t191 = -t87 / 0.2e1;
t189 = -t99 / 0.2e1;
t183 = pkin(4) * t105;
t177 = t93 * mrSges(4,2);
t175 = -qJ(5) - pkin(7);
t54 = -mrSges(6,2) * t99 + mrSges(6,3) * t86;
t55 = -mrSges(5,2) * t99 + mrSges(5,3) * t86;
t174 = t54 + t55;
t56 = mrSges(6,1) * t99 - mrSges(6,3) * t87;
t57 = mrSges(5,1) * t99 - mrSges(5,3) * t87;
t173 = -t56 - t57;
t89 = t133 * qJD(1);
t32 = t105 * t89 + t107 * t65;
t172 = t105 * t90 + t81 * t148;
t171 = t106 * t100 * t152 + t107 * t90;
t156 = t107 * t108;
t88 = t100 * t156;
t39 = t105 * t81 + t88;
t140 = mrSges(4,3) * t155;
t170 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t86 - mrSges(5,2) * t87 - t140;
t159 = qJ(5) * t106;
t62 = qJD(3) * t66;
t158 = t100 * t105;
t157 = t106 * t107;
t147 = qJD(5) * t107;
t139 = mrSges(4,3) * t154;
t16 = -t52 * mrSges(6,1) + t51 * mrSges(6,2);
t31 = -t105 * t65 + t107 * t89;
t134 = qJD(4) * t175;
t1 = pkin(4) * t135 - qJ(5) * t51 - qJD(5) * t87 + t4;
t2 = qJ(5) * t52 + qJD(5) * t86 + t3;
t132 = -t1 * t105 + t107 * t2;
t129 = mrSges(5,1) * t107 - mrSges(5,2) * t105;
t127 = mrSges(6,1) * t107 - mrSges(6,2) * t105;
t114 = pkin(4) * t106 - qJ(5) * t156;
t112 = -t4 * mrSges(5,1) - t1 * mrSges(6,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t101 = -pkin(4) * t107 - pkin(3);
t98 = Ifges(5,3) * t135;
t97 = Ifges(6,3) * t135;
t96 = t175 * t107;
t95 = t175 * t105;
t94 = -qJD(3) * mrSges(4,2) + t139;
t73 = (t100 + t183) * t106;
t72 = -qJD(5) * t105 + t107 * t134;
t71 = t105 * t134 + t147;
t70 = t107 * t81;
t48 = Ifges(5,5) * t51;
t47 = Ifges(6,5) * t51;
t46 = Ifges(5,6) * t52;
t45 = Ifges(6,6) * t52;
t43 = t146 + (qJD(1) * t183 + t91) * t108;
t42 = pkin(4) * t113 + t100 * t150;
t40 = -mrSges(6,1) * t86 + mrSges(6,2) * t87;
t38 = -t108 * t158 + t70;
t37 = -mrSges(5,2) * t135 + mrSges(5,3) * t52;
t36 = -mrSges(6,2) * t135 + mrSges(6,3) * t52;
t35 = mrSges(5,1) * t135 - mrSges(5,3) * t51;
t34 = mrSges(6,1) * t135 - mrSges(6,3) * t51;
t30 = -t105 * t159 + t39;
t23 = -qJ(5) * t157 + t70 + (-pkin(4) - t158) * t108;
t22 = -qJ(5) * t105 * t154 + t32;
t21 = -t52 * pkin(4) + t62;
t18 = qJD(1) * t114 + t31;
t17 = -mrSges(5,1) * t52 + mrSges(5,2) * t51;
t11 = -qJD(4) * t39 + t171;
t10 = (-t106 * t151 - t108 * t149) * t100 + t172;
t6 = (-qJ(5) * qJD(4) - qJD(3) * t100) * t157 + (-qJD(5) * t106 + (-qJ(5) * qJD(3) - qJD(4) * t100) * t108) * t105 + t172;
t5 = -t106 * t147 + t114 * qJD(3) + (-t88 + (-t81 + t159) * t105) * qJD(4) + t171;
t12 = [t10 * t55 + t11 * t57 + t73 * t16 + t23 * t34 + t30 * t36 + t38 * t35 + t39 * t37 + t42 * t40 + t5 * t56 + t6 * t54 + m(5) * (t20 * t10 + t19 * t11 + t3 * t39 + t4 * t38) + m(6) * (t1 * t23 + t2 * t30 + t21 * t73 + t33 * t42 + t5 * t7 + t6 * t9) + (-t97 / 0.2e1 - t98 / 0.2e1 - t47 / 0.2e1 - t48 / 0.2e1 - t45 / 0.2e1 - t46 / 0.2e1 + t61 * mrSges(4,3) - t142 * t52 - t143 * t51 + (0.3e1 / 0.2e1 * t102 + 0.2e1 * t177 + t137 + t210) * qJD(3) + t112) * t108 + (t21 * t126 + (-t1 * t107 - t105 * t2) * mrSges(6,3) + (-t105 * t3 - t107 * t4) * mrSges(5,3) + (mrSges(4,3) + t128) * t62 + (t136 - t196) * qJD(3) + (t59 * t129 + t33 * t127 + (t105 * t7 - t107 * t9) * mrSges(6,3) + (t105 * t19 - t107 * t20) * mrSges(5,3) + t203 * t193 + t202 * t191 + t204 * t189 - t205 * t107 / 0.2e1) * qJD(4) + (t138 * mrSges(4,1) + (-0.3e1 / 0.2e1 * Ifges(4,4) + t164 / 0.2e1 - t162 / 0.2e1 + t165 / 0.2e1 - t163 / 0.2e1) * t106 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t141) * t108) * t145 + t212 * t195 + t213 * t194 + t215 * t184) * t106 + (qJD(4) * t199 + t206) * t106 * t187 + ((m(4) * t61 + (-m(4) * t65 + m(5) * t59 - t170) * qJD(3)) * t108 + (t17 + (m(4) + m(5)) * t62 + (-m(4) * t66 - t94) * qJD(3)) * t106) * t100; (-t16 - t17 + (t105 * t173 + t107 * t174 - t139 + t94) * qJD(3) + m(5) * (t151 * t20 - t152 * t19 - t62) + m(6) * (t151 * t9 - t152 * t7 - t21)) * t108 + ((t36 + t37) * t107 + (-t34 - t35) * t105 + (-t105 * t174 + t107 * t173) * qJD(4) + (t40 - t140 - t170) * qJD(3) + m(5) * (qJD(3) * t59 - t148 * t19 - t149 * t20 + t131) + m(6) * (qJD(3) * t33 - t148 * t7 - t149 * t9 + t132)) * t106; m(5) * (-pkin(3) * t62 + pkin(7) * t131) + t202 * t195 + t203 * t194 + t215 * t105 / 0.2e1 + ((Ifges(4,4) * t209 + t204 * t216 + t136 + t196) * t106 + (-t177 - t102 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t155 + t137 - t210) * t108) * qJD(1) + ((m(6) * t33 + t40) * t183 + (-m(5) * t115 - t105 * t55 - t107 * t57) * pkin(7) - t197) * qJD(4) + t101 * t16 - t65 * t94 + t95 * t34 - t96 * t36 - m(5) * (t19 * t31 + t20 * t32 + t59 * t66) - m(6) * (t18 * t7 + t22 * t9 + t33 * t43) + (t72 - t18) * t56 - t61 * mrSges(4,2) - t32 * t55 - t31 * t57 - t43 * t40 - pkin(3) * t17 + (t71 - t22) * t54 + t206 * t184 + t132 * mrSges(6,3) + (-t105 * t35 + t107 * t37) * pkin(7) + t131 * mrSges(5,3) - t21 * t127 + (-mrSges(4,1) - t129) * t62 + t170 * t66 + m(6) * (t1 * t95 + t101 * t21 - t2 * t96 + t7 * t72 + t71 * t9); (-t40 * t87 + t34) * pkin(4) + (t7 * t86 + t87 * t9) * mrSges(6,3) + (t19 * t86 + t20 * t87) * mrSges(5,3) - t112 - t59 * (mrSges(5,1) * t87 + mrSges(5,2) * t86) - t33 * (mrSges(6,1) * t87 + mrSges(6,2) * t86) + t97 + t98 + t47 + t48 + t45 + t46 - t8 * t54 - t19 * t55 + t9 * t56 + t20 * t57 + (-(-t7 + t8) * t9 + (-t33 * t87 + t1) * pkin(4)) * m(6) + (t220 * t86 - t211) * t191 + t205 * t190 + (-t207 * t87 + t208 * t86) * t189 + (-t219 * t87 + t199 + t214) * t193; -t86 * t54 + t87 * t56 + 0.2e1 * (t21 / 0.2e1 + t7 * t190 + t9 * t193) * m(6) + t16;];
tauc = t12(:);

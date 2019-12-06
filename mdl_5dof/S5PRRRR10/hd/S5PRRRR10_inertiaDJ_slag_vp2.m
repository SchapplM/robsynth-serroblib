% Calculate time derivative of joint inertia matrix for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:19
% EndTime: 2019-12-05 17:23:27
% DurationCPUTime: 2.90s
% Computational Cost: add. (2901->406), mult. (8642->635), div. (0->0), fcn. (8168->12), ass. (0->177)
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t130 = sin(qJ(4));
t164 = qJD(5) * t130;
t134 = cos(qJ(4));
t166 = qJD(4) * t134;
t139 = -t129 * t164 + t133 * t166;
t211 = t129 / 0.2e1;
t194 = t133 / 0.2e1;
t125 = sin(pkin(6));
t131 = sin(qJ(3));
t176 = t125 * t131;
t118 = pkin(8) * t176;
t127 = cos(pkin(6));
t135 = cos(qJ(3));
t192 = pkin(2) * t135;
t89 = t127 * t192 - t118;
t136 = cos(qJ(2));
t171 = t135 * t136;
t132 = sin(qJ(2));
t174 = t131 * t132;
t210 = t127 * t171 - t174;
t209 = -m(5) * pkin(3) - mrSges(5,1) * t134 + mrSges(5,2) * t130 - mrSges(4,1);
t175 = t125 * t135;
t90 = t127 * t131 * pkin(2) + pkin(8) * t175;
t77 = pkin(9) * t127 + t90;
t78 = (-pkin(3) * t135 - pkin(9) * t131 - pkin(2)) * t125;
t185 = t130 * t78 + t134 * t77;
t169 = qJD(3) * t125;
t83 = (pkin(3) * t131 - pkin(9) * t135) * t169;
t84 = t89 * qJD(3);
t20 = -qJD(4) * t185 - t130 * t84 + t134 * t83;
t107 = -mrSges(6,1) * t133 + mrSges(6,2) * t129;
t208 = -m(6) * pkin(4) - mrSges(5,1) + t107;
t207 = 0.2e1 * m(6);
t206 = 0.2e1 * pkin(9);
t205 = -2 * mrSges(4,3);
t88 = t127 * t130 + t134 * t176;
t141 = t129 * t175 - t133 * t88;
t168 = qJD(3) * t131;
t156 = t125 * t168;
t155 = t135 * t169;
t87 = -t134 * t127 + t130 * t176;
t60 = -qJD(4) * t87 + t134 * t155;
t28 = qJD(5) * t141 - t129 * t60 + t133 * t156;
t203 = t28 / 0.2e1;
t62 = -t129 * t88 - t133 * t175;
t202 = t62 / 0.2e1;
t201 = -t141 / 0.2e1;
t181 = Ifges(6,4) * t129;
t147 = Ifges(6,1) * t133 - t181;
t81 = -Ifges(6,5) * t134 + t130 * t147;
t200 = t81 / 0.2e1;
t126 = sin(pkin(5));
t128 = cos(pkin(5));
t172 = t132 * t135;
t173 = t131 * t136;
t140 = t127 * t173 + t172;
t56 = t126 * t140 + t128 * t176;
t86 = -t125 * t126 * t136 + t128 * t127;
t35 = t130 * t56 - t86 * t134;
t170 = qJD(2) * t126;
t157 = t132 * t170;
t150 = t125 * t157;
t34 = t128 * t155 + (t210 * qJD(3) + (-t127 * t174 + t171) * qJD(2)) * t126;
t36 = t130 * t86 + t134 * t56;
t9 = qJD(4) * t36 + t130 * t34 - t134 * t150;
t199 = t35 * t9;
t198 = Ifges(6,5) * t211 + Ifges(6,6) * t194;
t180 = Ifges(6,4) * t133;
t112 = Ifges(6,1) * t129 + t180;
t197 = t112 / 0.2e1;
t196 = -t129 / 0.2e1;
t195 = -t133 / 0.2e1;
t48 = mrSges(5,1) * t156 - mrSges(5,3) * t60;
t27 = qJD(5) * t62 + t129 * t156 + t133 * t60;
t8 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t193 = t8 - t48;
t191 = pkin(9) * t134;
t33 = t128 * t156 + (t140 * qJD(3) + (t127 * t172 + t173) * qJD(2)) * t126;
t55 = -t210 * t126 - t128 * t175;
t190 = t33 * t55;
t85 = t90 * qJD(3);
t189 = t55 * t85;
t188 = t9 * t130;
t32 = -mrSges(6,1) * t62 - mrSges(6,2) * t141;
t65 = -mrSges(5,1) * t175 - mrSges(5,3) * t88;
t187 = t32 - t65;
t186 = -mrSges(4,1) * t127 + mrSges(5,1) * t87 + mrSges(5,2) * t88 + mrSges(4,3) * t176;
t184 = mrSges(6,3) * t130;
t183 = Ifges(5,4) * t130;
t182 = Ifges(5,4) * t134;
t179 = Ifges(6,6) * t129;
t10 = -qJD(4) * t35 + t130 * t150 + t134 * t34;
t178 = t10 * t134;
t167 = qJD(4) * t130;
t165 = qJD(5) * t129;
t163 = qJD(5) * t133;
t162 = qJD(5) * t134;
t61 = qJD(4) * t88 + t130 * t155;
t5 = Ifges(6,5) * t27 + Ifges(6,6) * t28 + Ifges(6,3) * t61;
t160 = Ifges(5,5) * t60 - Ifges(5,6) * t61 + Ifges(5,3) * t156;
t159 = Ifges(5,6) * t175;
t104 = (pkin(4) * t130 - pkin(10) * t134) * qJD(4);
t106 = -pkin(4) * t134 - pkin(10) * t130 - pkin(3);
t42 = t106 * t163 + t104 * t129 + (-t129 * t162 - t133 * t167) * pkin(9);
t73 = t106 * t133 - t129 * t191;
t152 = -qJD(5) * t73 + t42;
t43 = -t106 * t165 + t104 * t133 + (t129 * t167 - t133 * t162) * pkin(9);
t74 = t106 * t129 + t133 * t191;
t151 = -qJD(5) * t74 - t43;
t76 = t118 + (-pkin(3) - t192) * t127;
t37 = pkin(4) * t87 - pkin(10) * t88 + t76;
t39 = -pkin(10) * t175 + t185;
t11 = -t129 * t39 + t133 * t37;
t19 = t130 * t83 + t134 * t84 + t78 * t166 - t167 * t77;
t15 = pkin(10) * t156 + t19;
t24 = pkin(4) * t61 - pkin(10) * t60 + t85;
t1 = qJD(5) * t11 + t129 * t24 + t133 * t15;
t12 = t129 * t37 + t133 * t39;
t2 = -qJD(5) * t12 - t129 * t15 + t133 * t24;
t149 = t1 * t133 - t2 * t129;
t148 = mrSges(6,1) * t129 + mrSges(6,2) * t133;
t146 = -Ifges(6,2) * t129 + t180;
t110 = Ifges(6,2) * t133 + t181;
t18 = t129 * t55 + t133 * t36;
t17 = -t129 * t36 + t133 * t55;
t44 = -t130 * t77 + t134 * t78;
t143 = t166 * t35 + t188;
t22 = -Ifges(6,4) * t141 + Ifges(6,2) * t62 + Ifges(6,6) * t87;
t23 = -Ifges(6,1) * t141 + Ifges(6,4) * t62 + Ifges(6,5) * t87;
t142 = t194 * t23 + t196 * t22;
t138 = t129 * t166 + t130 * t163;
t50 = t139 * Ifges(6,5) - Ifges(6,6) * t138 + Ifges(6,3) * t167;
t124 = Ifges(5,5) * t166;
t123 = Ifges(6,5) * t163;
t116 = Ifges(4,5) * t155;
t113 = Ifges(5,1) * t130 + t182;
t111 = Ifges(5,2) * t134 + t183;
t103 = -mrSges(6,1) * t134 - t133 * t184;
t102 = mrSges(6,2) * t134 - t129 * t184;
t101 = (Ifges(5,1) * t134 - t183) * qJD(4);
t100 = t147 * qJD(5);
t99 = (-Ifges(5,2) * t130 + t182) * qJD(4);
t98 = t146 * qJD(5);
t97 = -Ifges(6,6) * t165 + t123;
t96 = (mrSges(5,1) * t130 + mrSges(5,2) * t134) * qJD(4);
t95 = t148 * qJD(5);
t94 = -mrSges(4,2) * t127 + mrSges(4,3) * t175;
t91 = t148 * t130;
t82 = (mrSges(4,1) * t131 + mrSges(4,2) * t135) * t169;
t80 = -Ifges(6,6) * t134 + t130 * t146;
t79 = -Ifges(6,3) * t134 + (Ifges(6,5) * t133 - t179) * t130;
t68 = -mrSges(6,2) * t167 - mrSges(6,3) * t138;
t67 = mrSges(6,1) * t167 - mrSges(6,3) * t139;
t64 = mrSges(5,2) * t175 - mrSges(5,3) * t87;
t54 = mrSges(6,1) * t138 + mrSges(6,2) * t139;
t52 = -t112 * t164 + (Ifges(6,5) * t130 + t134 * t147) * qJD(4);
t51 = -t110 * t164 + (Ifges(6,6) * t130 + t134 * t146) * qJD(4);
t49 = -mrSges(5,2) * t156 - mrSges(5,3) * t61;
t47 = Ifges(5,1) * t88 - Ifges(5,4) * t87 - Ifges(5,5) * t175;
t46 = Ifges(5,4) * t88 - Ifges(5,2) * t87 - t159;
t41 = mrSges(6,1) * t87 + mrSges(6,3) * t141;
t40 = -mrSges(6,2) * t87 + mrSges(6,3) * t62;
t38 = pkin(4) * t175 - t44;
t31 = mrSges(5,1) * t61 + mrSges(5,2) * t60;
t30 = Ifges(5,1) * t60 - Ifges(5,4) * t61 + Ifges(5,5) * t156;
t29 = Ifges(5,4) * t60 - Ifges(5,2) * t61 + Ifges(5,6) * t156;
t21 = -Ifges(6,5) * t141 + Ifges(6,6) * t62 + Ifges(6,3) * t87;
t16 = -pkin(4) * t156 - t20;
t14 = -mrSges(6,2) * t61 + mrSges(6,3) * t28;
t13 = mrSges(6,1) * t61 - mrSges(6,3) * t27;
t7 = Ifges(6,1) * t27 + Ifges(6,4) * t28 + Ifges(6,5) * t61;
t6 = Ifges(6,4) * t27 + Ifges(6,2) * t28 + Ifges(6,6) * t61;
t4 = qJD(5) * t17 + t10 * t133 + t129 * t33;
t3 = -qJD(5) * t18 - t10 * t129 + t133 * t33;
t25 = [0.2e1 * m(6) * (t17 * t3 + t18 * t4 + t199) + 0.2e1 * m(5) * (t10 * t36 + t190 + t199) + 0.2e1 * m(4) * (t150 * t86 + t34 * t56 + t190); t10 * t64 + t17 * t13 + t18 * t14 + t3 * t41 + t55 * t31 + t34 * t94 + t36 * t49 + t4 * t40 + t86 * t82 + t187 * t9 + t193 * t35 + t186 * t33 + (-mrSges(3,1) * t132 - mrSges(3,2) * t136) * t170 + ((-mrSges(4,1) * t135 + mrSges(4,2) * t131) * t150 + (-t131 * t56 + t135 * t55) * qJD(3) * mrSges(4,3)) * t125 + m(6) * (t1 * t18 + t11 * t3 + t12 * t4 + t16 * t35 + t17 * t2 + t38 * t9) + m(5) * (t10 * t185 + t19 * t36 - t20 * t35 + t33 * t76 - t44 * t9 + t189) + m(4) * (-pkin(2) * t125 ^ 2 * t157 - t33 * t89 + t34 * t90 + t56 * t84 + t189); 0.2e1 * m(4) * (t84 * t90 - t85 * t89) + 0.2e1 * m(5) * (t185 * t19 + t20 * t44 + t76 * t85) + 0.2e1 * t185 * t49 + (t5 - t29) * t87 + (t21 - t46) * t61 + 0.2e1 * t186 * t85 + (t1 * t12 + t11 * t2 + t16 * t38) * t207 + t127 * t116 - t141 * t7 + 0.2e1 * t11 * t13 + 0.2e1 * t12 * t14 + t27 * t23 + t28 * t22 + 0.2e1 * t16 * t32 + 0.2e1 * t38 * t8 + 0.2e1 * t1 * t40 + 0.2e1 * t2 * t41 + 0.2e1 * t44 * t48 + t60 * t47 + t62 * t6 + 0.2e1 * t19 * t64 + 0.2e1 * t20 * t65 + 0.2e1 * t76 * t31 + t88 * t30 + 0.2e1 * t84 * t94 + (-0.2e1 * pkin(2) * t82 - t135 * t160 + ((0.2e1 * Ifges(4,4) * t175 + Ifges(4,5) * t127 + t205 * t89) * t135 + (-0.2e1 * Ifges(4,4) * t176 + t90 * t205 + Ifges(5,5) * t88 - 0.2e1 * Ifges(4,6) * t127 - Ifges(5,6) * t87 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3)) * t175) * t131) * qJD(3)) * t125; -t34 * mrSges(4,2) + t4 * t102 + t3 * t103 + t17 * t67 + t18 * t68 + t35 * t54 + t55 * t96 + t9 * t91 + m(6) * (t17 * t43 + t18 * t42 + t3 * t73 + t4 * t74) + (m(6) * t143 / 0.2e1 + m(5) * (-t167 * t36 + t143 + t178) / 0.2e1) * t206 + (t178 + t188 + (-t130 * t36 + t134 * t35) * qJD(4)) * mrSges(5,3) + t209 * t33; ((-t44 * mrSges(5,3) + t47 / 0.2e1 + t142) * t134 + (-t185 * mrSges(5,3) + t21 / 0.2e1 - t46 / 0.2e1 + t159 / 0.2e1) * t130) * qJD(4) + (t134 * t49 + t193 * t130 + (-t130 * t64 + t134 * t187) * qJD(4) + m(6) * (t130 * t16 + t166 * t38) + m(5) * (-t20 * t130 + t19 * t134 - t166 * t44 - t167 * t185)) * pkin(9) + t27 * t200 + t52 * t201 + t51 * t202 + t80 * t203 + t209 * t85 + (t19 * mrSges(5,3) - t5 / 0.2e1 + t29 / 0.2e1) * t134 + t116 + m(6) * (t1 * t74 + t11 * t43 + t12 * t42 + t2 * t73) + (t50 / 0.2e1 - t99 / 0.2e1) * t87 - pkin(3) * t31 + t42 * t40 + t43 * t41 + t38 * t54 + t11 * t67 + t12 * t68 + t73 * t13 + t74 * t14 - t84 * mrSges(4,2) + t16 * t91 + t76 * t96 + t88 * t101 / 0.2e1 + t1 * t102 + t2 * t103 + t60 * t113 / 0.2e1 + (-t135 * t124 / 0.2e1 + (Ifges(5,5) * t130 / 0.2e1 + Ifges(5,6) * t134 / 0.2e1 - Ifges(4,6)) * t168) * t125 + (-t20 * mrSges(5,3) + t6 * t196 + t7 * t194 + t30 / 0.2e1 + (t195 * t22 + t196 * t23) * qJD(5)) * t130 + (t79 / 0.2e1 - t111 / 0.2e1) * t61; (t42 * t74 + t43 * t73) * t207 + 0.2e1 * t42 * t102 + 0.2e1 * t74 * t68 + 0.2e1 * t43 * t103 + 0.2e1 * t73 * t67 - 0.2e1 * pkin(3) * t96 + (-t50 + t99 + (-t129 * t80 + t133 * t81 + t206 * t91 + t113) * qJD(4)) * t134 + (t54 * t206 - t129 * t51 + t133 * t52 + t101 + (-t129 * t81 - t133 * t80) * qJD(5) + (pkin(9) ^ 2 * t134 * t207 - t111 + t79) * qJD(4)) * t130; -t10 * mrSges(5,2) + t35 * t95 + (m(6) * pkin(10) + mrSges(6,3)) * (-t3 * t129 + t4 * t133 + (-t129 * t18 - t133 * t17) * qJD(5)) + t208 * t9; -t19 * mrSges(5,2) + t20 * mrSges(5,1) + t38 * t95 + t87 * t97 / 0.2e1 + t98 * t202 + t100 * t201 + t16 * t107 + t61 * t198 + t110 * t203 + t27 * t197 + t7 * t211 + t6 * t194 + t142 * qJD(5) + (-m(6) * t16 - t8) * pkin(4) + ((-t11 * t133 - t12 * t129) * qJD(5) + t149) * mrSges(6,3) + (-t41 * t163 - t40 * t165 - t129 * t13 + t133 * t14 + m(6) * (-t11 * t163 - t12 * t165 + t149)) * pkin(10) + t160; -pkin(4) * t54 + t124 + (-t97 / 0.2e1 + t208 * qJD(4) * pkin(9)) * t134 + (t166 * t197 + qJD(5) * t200 + t51 / 0.2e1 + t152 * mrSges(6,3) + (m(6) * t152 - qJD(5) * t103 + t68) * pkin(10)) * t133 + (-t110 * t166 / 0.2e1 - qJD(5) * t80 / 0.2e1 + t52 / 0.2e1 + t151 * mrSges(6,3) + (m(6) * t151 - qJD(5) * t102 - t67) * pkin(10)) * t129 + (pkin(9) * t95 + t98 * t196 + t100 * t194 + (t110 * t195 + t112 * t196) * qJD(5) + (pkin(9) * mrSges(5,2) - Ifges(5,6) + t198) * qJD(4)) * t130; -0.2e1 * pkin(4) * t95 + t100 * t129 + t133 * t98 + (-t110 * t129 + t112 * t133) * qJD(5); mrSges(6,1) * t3 - mrSges(6,2) * t4; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t5; mrSges(6,1) * t43 - mrSges(6,2) * t42 + t50; t123 + (pkin(10) * t107 - t179) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t25(1), t25(2), t25(4), t25(7), t25(11); t25(2), t25(3), t25(5), t25(8), t25(12); t25(4), t25(5), t25(6), t25(9), t25(13); t25(7), t25(8), t25(9), t25(10), t25(14); t25(11), t25(12), t25(13), t25(14), t25(15);];
Mq = res;

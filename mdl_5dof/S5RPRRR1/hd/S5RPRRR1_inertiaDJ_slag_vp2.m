% Calculate time derivative of joint inertia matrix for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:55
% EndTime: 2019-07-18 13:25:00
% DurationCPUTime: 1.56s
% Computational Cost: add. (704->258), mult. (2263->431), div. (0->0), fcn. (1728->6), ass. (0->136)
t77 = cos(qJ(4));
t120 = qJD(4) * t77;
t78 = cos(qJ(3));
t123 = qJD(3) * t78;
t74 = sin(qJ(4));
t75 = sin(qJ(3));
t80 = t75 * t120 + t74 * t123;
t137 = t74 * t75;
t76 = cos(qJ(5));
t134 = t76 * t78;
t135 = t75 * t77;
t73 = sin(qJ(5));
t36 = -t73 * t135 - t134;
t37 = t76 * t135 - t73 * t78;
t13 = Ifges(6,5) * t37 + Ifges(6,6) * t36 + Ifges(6,3) * t137;
t138 = Ifges(5,6) * t78;
t141 = Ifges(5,4) * t77;
t88 = -Ifges(5,2) * t74 + t141;
t31 = t75 * t88 - t138;
t163 = t13 - t31;
t114 = qJ(2) * qJD(3);
t162 = -t78 * qJD(2) + t75 * t114;
t118 = qJD(5) * t74;
t161 = t73 * t118 - t76 * t120;
t117 = qJD(5) * t76;
t160 = -t74 * t117 - t73 * t120;
t68 = t74 ^ 2;
t71 = t77 ^ 2;
t126 = t68 + t71;
t91 = mrSges(6,1) * t73 + mrSges(6,2) * t76;
t38 = t91 * t74;
t159 = t126 * mrSges(5,3) + t74 * t38 - mrSges(4,2);
t122 = qJD(4) * t74;
t96 = qJD(3) * t77 - qJD(5);
t158 = t78 * t122 + t75 * t96;
t157 = 2 * m(6);
t92 = mrSges(5,1) * t74 + mrSges(5,2) * t77;
t156 = 0.2e1 * t92 * t75;
t155 = 0.2e1 * qJ(2);
t121 = qJD(4) * t75;
t110 = t74 * t121;
t136 = t75 * t76;
t97 = -qJD(5) * t77 + qJD(3);
t11 = t97 * t136 + (-t78 * t96 + t110) * t73;
t154 = t11 / 0.2e1;
t140 = Ifges(6,4) * t73;
t89 = Ifges(6,1) * t76 - t140;
t32 = -Ifges(6,5) * t77 + t74 * t89;
t153 = t32 / 0.2e1;
t152 = t36 / 0.2e1;
t151 = t37 / 0.2e1;
t139 = Ifges(6,4) * t76;
t55 = Ifges(6,1) * t73 + t139;
t150 = t55 / 0.2e1;
t149 = -t73 / 0.2e1;
t148 = -t74 / 0.2e1;
t147 = -t76 / 0.2e1;
t146 = t76 / 0.2e1;
t145 = t77 / 0.2e1;
t111 = t77 * t123;
t12 = t96 * t134 + (-t122 * t76 + t73 * t97) * t75;
t124 = qJD(3) * t75;
t144 = mrSges(6,1) * t11 - mrSges(6,2) * t12 + mrSges(5,1) * t124 + (t110 - t111) * mrSges(5,3);
t143 = mrSges(6,3) * t74;
t142 = Ifges(5,4) * t74;
t47 = mrSges(5,2) * t78 - mrSges(5,3) * t137;
t133 = t77 * t47;
t132 = t77 * t78;
t131 = mrSges(6,1) * t76 - mrSges(6,2) * t73 + mrSges(5,1);
t130 = mrSges(5,1) * t78 - mrSges(6,1) * t36 + mrSges(6,2) * t37 + mrSges(5,3) * t135;
t93 = mrSges(5,1) * t77 - mrSges(5,2) * t74;
t129 = -t93 - mrSges(4,1);
t128 = -Ifges(5,5) * t111 - Ifges(5,3) * t124;
t127 = t73 ^ 2 + t76 ^ 2;
t125 = qJ(2) * t78;
t119 = qJD(5) * t73;
t115 = qJ(2) * qJD(2);
t52 = Ifges(6,5) * t73 + Ifges(6,6) * t76;
t113 = t52 / 0.2e1 - Ifges(5,6);
t112 = 0.2e1 * t74;
t109 = t74 * t120;
t85 = Ifges(6,5) * t76 - Ifges(6,6) * t73;
t29 = -Ifges(6,3) * t77 + t74 * t85;
t54 = Ifges(5,2) * t77 + t142;
t101 = -t54 / 0.2e1 + t29 / 0.2e1;
t72 = t78 ^ 2;
t99 = t72 * t115;
t1 = Ifges(6,5) * t12 + Ifges(6,6) * t11 + t80 * Ifges(6,3);
t95 = 0.2e1 * t130;
t83 = t76 * t132 + t73 * t75;
t84 = t97 * t78;
t5 = t83 * qJD(2) + (-t158 * t76 + t73 * t84) * qJ(2);
t82 = -t73 * t132 + t136;
t6 = t82 * qJD(2) + (t158 * t73 + t76 * t84) * qJ(2);
t94 = t5 * t76 - t6 * t73;
t90 = Ifges(5,1) * t77 - t142;
t56 = Ifges(5,1) * t74 + t141;
t87 = -Ifges(6,2) * t73 + t139;
t53 = Ifges(6,2) * t76 + t140;
t86 = Ifges(5,5) * t74 + Ifges(5,6) * t77;
t14 = Ifges(6,4) * t37 + Ifges(6,2) * t36 + Ifges(6,6) * t137;
t15 = Ifges(6,1) * t37 + Ifges(6,4) * t36 + Ifges(6,5) * t137;
t81 = t14 * t149 + t15 * t146;
t79 = qJ(2) ^ 2;
t69 = t75 ^ 2;
t66 = Ifges(5,5) * t120;
t64 = t69 * t115;
t57 = t68 * t99;
t48 = -mrSges(6,1) * t77 - t76 * t143;
t46 = mrSges(6,2) * t77 - t73 * t143;
t45 = t90 * qJD(4);
t44 = t89 * qJD(5);
t43 = t88 * qJD(4);
t42 = t87 * qJD(5);
t41 = t85 * qJD(5);
t40 = t91 * qJD(5);
t35 = t83 * qJ(2);
t34 = t82 * qJ(2);
t33 = -t78 * Ifges(5,5) + t75 * t90;
t30 = -Ifges(6,6) * t77 + t74 * t87;
t28 = -mrSges(5,2) * t124 - mrSges(5,3) * t80;
t26 = -mrSges(6,2) * t122 + t160 * mrSges(6,3);
t25 = mrSges(6,1) * t122 + t161 * mrSges(6,3);
t24 = mrSges(6,1) * t137 - mrSges(6,3) * t37;
t23 = -mrSges(6,2) * t137 + mrSges(6,3) * t36;
t22 = t160 * mrSges(6,1) + t161 * mrSges(6,2);
t20 = -t56 * t121 + (Ifges(5,5) * t75 + t78 * t90) * qJD(3);
t19 = -t55 * t118 + (Ifges(6,5) * t74 + t77 * t89) * qJD(4);
t18 = -t54 * t121 + (Ifges(5,6) * t75 + t78 * t88) * qJD(3);
t17 = -t53 * t118 + (Ifges(6,6) * t74 + t77 * t87) * qJD(4);
t16 = -t52 * t118 + (Ifges(6,3) * t74 + t77 * t85) * qJD(4);
t8 = mrSges(6,1) * t80 - mrSges(6,3) * t12;
t7 = -mrSges(6,2) * t80 + mrSges(6,3) * t11;
t3 = Ifges(6,1) * t12 + Ifges(6,4) * t11 + Ifges(6,5) * t80;
t2 = Ifges(6,4) * t12 + Ifges(6,2) * t11 + Ifges(6,6) * t80;
t4 = [t11 * t14 + t12 * t15 + t36 * t2 + 0.2e1 * t5 * t23 + 0.2e1 * t6 * t24 + t37 * t3 + 0.2e1 * t34 * t8 + 0.2e1 * t35 * t7 + (m(3) * t155 + (2 * mrSges(3,3)) + 0.2e1 * (t69 + t72) * mrSges(4,3)) * qJD(2) + (t109 * t72 * t79 + t34 * t6 + t35 * t5 + t57) * t157 + 0.2e1 * m(5) * (t71 * t99 + t57 + t64) + 0.2e1 * m(4) * (t64 + t99) + ((t74 * t95 + 0.2e1 * t133) * qJD(2) + (0.2e1 * Ifges(4,4) * t78 + t77 * t33 + (t138 + t163) * t74) * qJD(3) + (qJD(3) * t156 + 0.2e1 * t77 * t28 - t144 * t112 + (-0.2e1 * t47 * t74 + t77 * t95) * qJD(4)) * qJ(2) + t128) * t78 + (qJD(2) * t156 + t77 * t20 + (t1 - t18) * t74 + (t75 * t93 * t155 + t163 * t77 - t74 * t33 + t78 * t86) * qJD(4) + ((Ifges(5,5) * t77 - Ifges(5,6) * t74 - 0.2e1 * Ifges(4,4)) * t75 + (-t130 * t112 - 0.2e1 * t133) * qJ(2) + (t92 * t155 - Ifges(5,3) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + 0.4e1 * (-m(6) * t68 / 0.2e1 + m(5) * (0.1e1 - t126) / 0.2e1) * t79) * t78) * qJD(3)) * t75; t144 * t77 + (t75 * mrSges(4,1) + t78 * mrSges(4,2)) * qJD(3) + (m(6) * (t68 - t71) * t125 + (m(6) * (-t34 * t73 + t35 * t76) + t76 * t23 - t73 * t24 + t47) * t77) * qJD(4) + (m(6) * (-t117 * t34 - t119 * t35 + t162 * t77 + t94) - t23 * t119 + t76 * t7 - t24 * t117 - t73 * t8 + t28 + t130 * qJD(4)) * t74; (-0.1e1 + t127) * t109 * t157; t17 * t152 + t19 * t151 + t5 * t46 + t6 * t48 + t30 * t154 + t12 * t153 + t34 * t25 + t35 * t26 + (t18 / 0.2e1 - t1 / 0.2e1 + (t33 / 0.2e1 + t81) * qJD(4)) * t77 + (t3 * t146 + t2 * t149 + t20 / 0.2e1 + (t14 * t147 + t15 * t149) * qJD(5) + (-t31 / 0.2e1 + t13 / 0.2e1) * qJD(4)) * t74 + (Ifges(5,6) * t122 / 0.2e1 - t66 / 0.2e1 + (t120 * t38 - t74 * t22) * qJ(2) + t159 * qJD(2) + (t129 * qJ(2) + t101 * t74 + t56 * t145 + Ifges(4,5)) * qJD(3)) * t78 + (t45 * t145 + t43 * t148 + t74 * t16 / 0.2e1 + t129 * qJD(2) + (qJ(2) * t92 + t101 * t77 + t56 * t148) * qJD(4) + (t86 / 0.2e1 - Ifges(4,6) - t159 * qJ(2)) * qJD(3)) * t75; (t22 + (t46 * t76 - t48 * t73) * qJD(4)) * t77 + (qJD(4) * t38 - t25 * t73 + t26 * t76 + (-t46 * t73 - t48 * t76) * qJD(5)) * t74; (-t16 + t43 + (-t73 * t30 + t76 * t32 + t56) * qJD(4)) * t77 + (-t17 * t73 + t19 * t76 + t45 + (-t30 * t76 - t32 * t73) * qJD(5) + (t29 - t54) * qJD(4)) * t74; t2 * t146 + t73 * t3 / 0.2e1 + t42 * t152 + t44 * t151 + t53 * t154 + t12 * t150 + t81 * qJD(5) + ((-t34 * t76 - t35 * t73) * qJD(5) + t94) * mrSges(6,3) + (t162 * mrSges(5,2) + (t113 * t75 - t131 * t125) * qJD(4)) * t77 + ((-Ifges(5,5) * qJD(4) + t41 / 0.2e1 + t131 * t114) * t75 + (t113 * qJD(3) - t131 * qJD(2) + (mrSges(5,2) * qJD(4) + t40) * qJ(2)) * t78) * t74 - t128; -t40 * t77 + (-t131 * t74 + (t127 * mrSges(6,3) - mrSges(5,2)) * t77) * qJD(4); -t77 * t41 / 0.2e1 + t66 + (t120 * t150 + qJD(5) * t153 + t17 / 0.2e1) * t76 + (-t53 * t120 / 0.2e1 + t19 / 0.2e1 - qJD(5) * t30 / 0.2e1) * t73 + (t44 * t146 + t42 * t149 + (t53 * t147 + t55 * t149) * qJD(5) + t113 * qJD(4)) * t74; t42 * t76 + t44 * t73 + (-t53 * t73 + t55 * t76) * qJD(5); mrSges(6,1) * t6 - t5 * mrSges(6,2) + t1; t22; t16; t41; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq  = res;

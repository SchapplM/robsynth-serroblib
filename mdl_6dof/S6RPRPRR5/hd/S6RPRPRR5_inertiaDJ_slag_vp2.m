% Calculate time derivative of joint inertia matrix for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:48:01
% EndTime: 2019-03-09 03:48:06
% DurationCPUTime: 2.37s
% Computational Cost: add. (4387->299), mult. (9360->424), div. (0->0), fcn. (9383->8), ass. (0->132)
t132 = (mrSges(5,2) + mrSges(4,3));
t173 = 2 * t132;
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t129 = t84 ^ 2 + t87 ^ 2;
t116 = t129 * mrSges(7,3);
t88 = cos(qJ(5));
t172 = (mrSges(6,2) - t116) * t88;
t149 = cos(qJ(3));
t111 = t149 * qJD(2);
t131 = pkin(7) + qJ(2);
t83 = cos(pkin(10));
t86 = sin(qJ(3));
t137 = t86 * t83;
t64 = t131 * t83;
t53 = t149 * t64;
t82 = sin(pkin(10));
t32 = (-qJD(3) * t131 * t86 + t111) * t82 + qJD(2) * t137 + qJD(3) * t53;
t95 = t149 * t83 - t86 * t82;
t46 = t95 * qJD(3);
t170 = -t46 * pkin(8) + t32;
t65 = -t87 * mrSges(7,1) + mrSges(7,2) * t84;
t133 = mrSges(6,1) - t65;
t161 = -m(7) * pkin(5) - t133;
t150 = -t83 * pkin(2) - pkin(1);
t168 = 0.2e1 * t150;
t166 = m(7) * pkin(9);
t85 = sin(qJ(5));
t89 = -pkin(3) - pkin(4);
t62 = -t85 * qJ(4) + t88 * t89;
t41 = t88 * qJD(4) + qJD(5) * t62;
t165 = t41 * mrSges(6,2);
t63 = t88 * qJ(4) + t85 * t89;
t164 = t129 * t88;
t55 = t149 * t82 + t137;
t102 = -t55 * t85 - t88 * t95;
t35 = -pkin(3) * t95 - t55 * qJ(4) + t150;
t30 = pkin(4) * t95 - t35;
t37 = t55 * t88 - t85 * t95;
t13 = -pkin(5) * t102 - pkin(9) * t37 + t30;
t117 = t131 * t82;
t39 = -t86 * t117 + t53;
t34 = -pkin(8) * t95 + t39;
t52 = t149 * t117;
t38 = t64 * t86 + t52;
t98 = -pkin(8) * t55 + t38;
t17 = t88 * t34 + t85 * t98;
t10 = t13 * t84 + t17 * t87;
t9 = t13 * t87 - t17 * t84;
t163 = t10 * t87 - t84 * t9;
t120 = qJD(6) * t87;
t16 = t85 * t34 - t88 * t98;
t47 = t55 * qJD(3);
t21 = qJD(5) * t102 + t46 * t88 + t47 * t85;
t125 = qJD(5) * t37;
t22 = t46 * t85 - t88 * t47 + t125;
t104 = mrSges(7,1) * t84 + mrSges(7,2) * t87;
t56 = t104 * qJD(6);
t127 = qJD(5) * t16;
t31 = -qJD(3) * t52 + t83 * t111 + (-qJD(2) * t82 - qJD(3) * t64) * t86;
t28 = pkin(8) * t47 + t31;
t6 = t170 * t85 + t88 * t28 - t127;
t66 = Ifges(7,5) * t84 + Ifges(7,6) * t87;
t121 = qJD(6) * t84;
t71 = Ifges(7,6) * t121;
t162 = (t66 / 0.2e1 - Ifges(6,6)) * t22 + Ifges(6,5) * t21 + t16 * t56 - t102 * (Ifges(7,5) * t120 - t71) / 0.2e1 - t6 * mrSges(6,2);
t58 = Ifges(7,4) * t120 - Ifges(7,2) * t121;
t59 = Ifges(7,1) * t120 - Ifges(7,4) * t121;
t147 = Ifges(7,4) * t84;
t67 = Ifges(7,2) * t87 + t147;
t146 = Ifges(7,4) * t87;
t68 = Ifges(7,1) * t84 + t146;
t90 = -(t84 * t67 - t87 * t68) * qJD(6) + t87 * t58 + t84 * t59;
t160 = 2 * m(6);
t159 = 0.2e1 * m(7);
t126 = qJD(5) * t17;
t7 = -t88 * t170 + t85 * t28 + t126;
t158 = 0.2e1 * t7;
t99 = qJ(4) * t46 + qJD(4) * t55;
t26 = t47 * t89 + t99;
t157 = 0.2e1 * t26;
t156 = -0.2e1 * t56;
t155 = -t37 / 0.2e1;
t154 = -t67 / 0.2e1;
t153 = t16 * t7;
t148 = mrSges(7,3) * t37;
t145 = Ifges(7,5) * t87;
t42 = t85 * qJD(4) + qJD(5) * t63;
t143 = t16 * t42;
t142 = t41 * t85;
t141 = t42 * t88;
t138 = t84 * t21;
t136 = t87 * t21;
t135 = t88 * t56;
t134 = qJD(6) / 0.2e1;
t130 = Ifges(7,5) * t136 + Ifges(7,3) * t22;
t128 = qJD(6) * t9;
t124 = qJD(5) * t88;
t123 = qJD(6) * t10;
t61 = -pkin(9) + t63;
t122 = qJD(6) * t61;
t118 = t37 * t121;
t115 = t129 * t41;
t114 = -Ifges(7,6) * t84 - (2 * Ifges(6,4));
t8 = pkin(5) * t22 - pkin(9) * t21 + t26;
t1 = t6 * t87 + t8 * t84 + t128;
t113 = -t1 + t128;
t112 = t39 * t31 + t32 * t38;
t2 = -t6 * t84 + t8 * t87 - t123;
t110 = t2 + t123;
t108 = -2 * Ifges(4,4) + 2 * Ifges(5,5);
t96 = t118 - t136;
t97 = t120 * t37 + t138;
t3 = -Ifges(7,4) * t96 - Ifges(7,2) * t97 + t22 * Ifges(7,6);
t107 = t3 / 0.2e1 + t21 * t68 / 0.2e1;
t4 = -Ifges(7,1) * t96 - Ifges(7,4) * t97 + t22 * Ifges(7,5);
t106 = t4 / 0.2e1 + t21 * t154;
t105 = t22 * mrSges(6,1) + t21 * mrSges(6,2);
t23 = mrSges(7,2) * t102 - t148 * t84;
t24 = -mrSges(7,1) * t102 - t148 * t87;
t103 = -t87 * t23 + t84 * t24;
t60 = pkin(5) - t62;
t44 = t47 * mrSges(5,1);
t43 = t46 * mrSges(4,2);
t29 = pkin(3) * t47 - t99;
t19 = t104 * t37;
t15 = -t102 * Ifges(7,5) + (Ifges(7,1) * t87 - t147) * t37;
t14 = -t102 * Ifges(7,6) + (-Ifges(7,2) * t84 + t146) * t37;
t12 = -mrSges(7,2) * t22 - mrSges(7,3) * t97;
t11 = mrSges(7,1) * t22 + mrSges(7,3) * t96;
t5 = mrSges(7,1) * t97 - mrSges(7,2) * t96;
t18 = [t19 * t158 + (t1 * t10 + t2 * t9 + t153) * t159 + (t17 * t6 + t26 * t30 + t153) * t160 + t15 * t136 + (mrSges(6,2) * t157 + mrSges(6,3) * t158 + 0.2e1 * Ifges(6,1) * t21 - t84 * t3 + t87 * t4 + (t114 + t145) * t22 + (t102 * t66 - t87 * t14 - t84 * t15) * qJD(6)) * t37 + 0.2e1 * t29 * (-mrSges(5,1) * t95 - mrSges(5,3) * t55) - (mrSges(6,1) * t157 - 0.2e1 * t6 * mrSges(6,3) + ((2 * Ifges(6,2)) + Ifges(7,3)) * t22 + t114 * t21 + t130) * t102 - t14 * t138 + 0.2e1 * m(5) * (t29 * t35 + t112) + 0.2e1 * m(4) * t112 + 0.2e1 * t30 * t105 + 0.2e1 * t35 * t44 + 0.2e1 * (t16 * t21 - t17 * t22) * mrSges(6,3) + (mrSges(4,1) * t168 + t55 * t108 - 0.2e1 * (Ifges(5,3) + Ifges(4,2)) * t95 - 0.2e1 * t132 * t39) * t47 + 0.2e1 * t1 * t23 + 0.2e1 * t2 * t24 + 0.2e1 * t16 * t5 + 0.2e1 * t9 * t11 + 0.2e1 * t10 * t12 + (-0.2e1 * t35 * mrSges(5,3) - t95 * t108 + 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t55 + t38 * t173) * t46 + t43 * t168 + (t31 * t95 + t32 * t55) * t173 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t82 ^ 2 + t83 ^ 2) * qJD(2); t47 * mrSges(4,1) - t46 * mrSges(5,3) - t87 * t11 - t84 * t12 + t43 + t44 + t103 * qJD(6) + m(7) * (-qJD(6) * t163 - t1 * t84 - t2 * t87) - m(6) * t26 + m(5) * t29 - t105; 0; t60 * t5 + t42 * t19 + t133 * t7 + (Ifges(5,6) - Ifges(4,6)) * t47 + (Ifges(5,4) + Ifges(4,5)) * t46 + (-mrSges(5,1) - mrSges(4,1)) * t32 + (-mrSges(4,2) + mrSges(5,3)) * t31 + (-pkin(3) * t46 - qJ(4) * t47 + qJD(4) * t95) * mrSges(5,2) + m(7) * (t60 * t7 + t143) + m(5) * (-pkin(3) * t32 + qJ(4) * t31 + qJD(4) * t39) + m(6) * (t17 * t41 + t6 * t63 - t62 * t7 + t143) + (t102 * t41 - t21 * t62 - t22 * t63 + t37 * t42) * mrSges(6,3) + (-t24 * t122 - qJD(6) * t15 / 0.2e1 + t41 * t23 + t61 * t12 + m(7) * (t1 * t61 + t10 * t41 - t122 * t9) + (t67 * t134 - t59 / 0.2e1) * t37 + t113 * mrSges(7,3) - t107) * t87 + (-t23 * t122 + t14 * t134 - t41 * t24 - t61 * t11 + m(7) * (-t10 * t122 - t2 * t61 - t41 * t9) + (t68 * t134 + t58 / 0.2e1) * t37 + t110 * mrSges(7,3) - t106) * t84 - t162; 0; t60 * t156 + 0.2e1 * t165 + (t115 * t61 + t42 * t60) * t159 + (t41 * t63 - t42 * t62) * t160 + t90 + 0.2e1 * t133 * t42 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4) - 0.2e1 * t116 * t41; t46 * mrSges(5,2) + m(5) * t32 + (-t21 * mrSges(6,3) - t5 - m(7) * t7 + m(6) * (-t7 + t126) + (m(7) * t163 + mrSges(6,3) * t102 - t103) * qJD(5)) * t88 + (qJD(5) * t19 - t84 * t11 + t87 * t12 + (-t84 * t23 - t87 * t24) * qJD(6) + (-t22 + t125) * mrSges(6,3) + m(7) * (t1 * t87 - t10 * t121 - t120 * t9 - t2 * t84 + t127) + m(6) * (t6 + t127)) * t85; 0; t135 + m(7) * (t129 * t142 - t141) + m(6) * (-t141 + t142) + (t133 * t85 + t172 + m(7) * (t164 * t61 + t60 * t85) + m(6) * (-t62 * t85 + t63 * t88)) * qJD(5); (-0.1e1 + t129) * t85 * t124 * t159; -pkin(5) * t5 + t161 * t7 + (t1 * mrSges(7,3) + t37 * t59 / 0.2e1 + (t37 * t154 - t9 * mrSges(7,3) + t15 / 0.2e1) * qJD(6) + (-m(7) * t113 - qJD(6) * t24 + t12) * pkin(9) + t107) * t87 + (-t2 * mrSges(7,3) + t58 * t155 + (t68 * t155 - t10 * mrSges(7,3) - t14 / 0.2e1) * qJD(6) + (-m(7) * t110 - qJD(6) * t23 - t11) * pkin(9) + t106) * t84 + t162; 0; -t165 + (pkin(5) + t60) * t56 + (mrSges(7,3) + t166) * t115 + t161 * t42 - t90; -t135 + (t161 * t85 + t164 * t166 - t172) * qJD(5); pkin(5) * t156 + t90; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,5) * t118 - Ifges(7,6) * t97 + t130; t56; t71 - t104 * t41 + (t61 * t65 - t145) * qJD(6); (t121 * t85 - t124 * t87) * mrSges(7,2) + (-t120 * t85 - t124 * t84) * mrSges(7,1); -t71 + (pkin(9) * t65 + t145) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;

% Calculate time derivative of joint inertia matrix for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:32
% EndTime: 2019-03-09 03:04:35
% DurationCPUTime: 1.96s
% Computational Cost: add. (2392->298), mult. (5150->418), div. (0->0), fcn. (4514->8), ass. (0->131)
t138 = Ifges(7,4) + Ifges(6,5);
t89 = sin(qJ(5));
t91 = cos(qJ(5));
t128 = t89 ^ 2 + t91 ^ 2;
t86 = sin(pkin(10));
t87 = cos(pkin(10));
t90 = sin(qJ(3));
t92 = cos(qJ(3));
t61 = t86 * t90 - t87 * t92;
t58 = t61 * qJD(3);
t167 = t128 * t58;
t166 = Ifges(7,2) + Ifges(6,3);
t125 = qJD(5) * t91;
t141 = t89 * t58;
t62 = t86 * t92 + t87 * t90;
t97 = t62 * t125 - t141;
t126 = qJD(5) * t89;
t139 = t91 * t58;
t96 = t62 * t126 + t139;
t115 = -cos(pkin(9)) * pkin(1) - pkin(2);
t69 = -pkin(3) * t92 + t115;
t34 = t61 * pkin(4) - t62 * pkin(8) + t69;
t78 = sin(pkin(9)) * pkin(1) + pkin(7);
t127 = qJ(4) + t78;
t113 = t127 * t90;
t60 = t127 * t92;
t36 = -t113 * t86 + t87 * t60;
t163 = t89 * t34 + t91 * t36;
t165 = qJD(5) * t163;
t111 = qJD(3) * t127;
t46 = qJD(4) * t92 - t111 * t90;
t94 = -t90 * qJD(4) - t111 * t92;
t24 = t87 * t46 + t86 * t94;
t120 = pkin(3) * qJD(3) * t90;
t57 = t62 * qJD(3);
t32 = pkin(4) * t57 + pkin(8) * t58 + t120;
t4 = -t24 * t89 + t32 * t91 - t165;
t2 = -pkin(5) * t57 - t4;
t164 = m(7) * t2;
t162 = Ifges(7,6) * t126 + t138 * t125;
t103 = pkin(5) * t91 + qJ(6) * t89;
t124 = qJD(6) * t91;
t161 = qJD(5) * t103 - t124;
t13 = t34 * t91 - t36 * t89;
t142 = t62 * t91;
t40 = mrSges(6,1) * t61 - mrSges(6,3) * t142;
t41 = -mrSges(7,1) * t61 + mrSges(7,2) * t142;
t132 = -t40 + t41;
t20 = -mrSges(6,2) * t57 - mrSges(6,3) * t97;
t21 = -mrSges(7,2) * t97 + mrSges(7,3) * t57;
t136 = t20 + t21;
t3 = t34 * t125 - t126 * t36 + t91 * t24 + t89 * t32;
t160 = t132 * qJD(5) + m(6) * (-t13 * qJD(5) + t3) + t136;
t143 = t62 * t89;
t39 = -mrSges(6,2) * t61 - mrSges(6,3) * t143;
t42 = -mrSges(7,2) * t143 + mrSges(7,3) * t61;
t133 = t39 + t42;
t18 = mrSges(6,1) * t57 + mrSges(6,3) * t96;
t19 = -t57 * mrSges(7,1) - t96 * mrSges(7,2);
t137 = t18 - t19;
t159 = -t133 * qJD(5) + m(6) * (-t4 - t165) - t137;
t158 = 2 * m(5);
t157 = -2 * mrSges(5,3);
t156 = -2 * Ifges(5,4);
t23 = t46 * t86 - t87 * t94;
t155 = 0.2e1 * t23;
t35 = t87 * t113 + t60 * t86;
t154 = 0.2e1 * t35;
t152 = Ifges(6,4) * t89;
t151 = Ifges(6,4) * t91;
t150 = Ifges(7,5) * t89;
t149 = Ifges(7,5) * t91;
t148 = t23 * t35;
t147 = t57 * t61;
t146 = t58 * t62;
t144 = t61 * Ifges(6,6);
t104 = Ifges(7,3) * t89 + t149;
t25 = Ifges(7,6) * t61 + t104 * t62;
t105 = -Ifges(6,2) * t89 + t151;
t26 = t105 * t62 + t144;
t135 = t25 - t26;
t106 = Ifges(7,1) * t91 + t150;
t27 = Ifges(7,4) * t61 + t106 * t62;
t107 = Ifges(6,1) * t91 - t152;
t28 = Ifges(6,5) * t61 + t107 * t62;
t134 = t27 + t28;
t77 = pkin(3) * t86 + pkin(8);
t131 = t77 * t167;
t71 = -t91 * mrSges(6,1) + t89 * mrSges(6,2);
t130 = t71 - mrSges(5,1);
t72 = -Ifges(7,3) * t91 + t150;
t73 = Ifges(6,2) * t91 + t152;
t117 = t72 / 0.2e1 - t73 / 0.2e1;
t74 = Ifges(7,1) * t89 - t149;
t75 = Ifges(6,1) * t89 + t151;
t116 = t74 / 0.2e1 + t75 / 0.2e1;
t79 = -pkin(3) * t87 - pkin(4);
t53 = t58 * mrSges(5,2);
t114 = t57 * mrSges(5,1) - t53;
t112 = 0.2e1 * t120;
t110 = mrSges(4,1) * t90 + mrSges(4,2) * t92;
t109 = t89 * mrSges(6,1) + t91 * mrSges(6,2);
t70 = -t91 * mrSges(7,1) - t89 * mrSges(7,3);
t108 = t89 * mrSges(7,1) - t91 * mrSges(7,3);
t102 = pkin(5) * t89 - qJ(6) * t91;
t101 = t23 * t61 + t35 * t57;
t98 = t97 * Ifges(7,6) - t138 * t139 + t166 * t57;
t95 = t132 * t89 + t133 * t91;
t93 = m(7) * t124 + (-m(7) * t103 + t70 + t71) * qJD(5);
t68 = t107 * qJD(5);
t67 = t106 * qJD(5);
t66 = t105 * qJD(5);
t65 = t104 * qJD(5);
t64 = t109 * qJD(5);
t63 = t108 * qJD(5);
t59 = -t103 + t79;
t56 = -pkin(5) * t126 + qJ(6) * t125 + qJD(6) * t89;
t38 = t109 * t62;
t37 = t108 * t62;
t17 = t102 * t62 + t35;
t16 = mrSges(6,1) * t97 - mrSges(6,2) * t96;
t15 = mrSges(7,1) * t97 + mrSges(7,3) * t96;
t11 = -Ifges(6,1) * t96 - Ifges(6,4) * t97 + t57 * Ifges(6,5);
t10 = -Ifges(7,1) * t96 + t57 * Ifges(7,4) + Ifges(7,5) * t97;
t9 = -Ifges(6,4) * t96 - Ifges(6,2) * t97 + t57 * Ifges(6,6);
t8 = -Ifges(7,5) * t96 + t57 * Ifges(7,6) + Ifges(7,3) * t97;
t7 = -pkin(5) * t61 - t13;
t6 = qJ(6) * t61 + t163;
t5 = -t102 * t58 + t161 * t62 + t23;
t1 = qJ(6) * t57 + qJD(6) * t61 + t3;
t12 = [t36 * t57 * t157 + t16 * t154 + 0.2e1 * t5 * t37 + t38 * t155 + 0.2e1 * t3 * t39 + 0.2e1 * t4 * t40 + 0.2e1 * t2 * t41 + 0.2e1 * t1 * t42 + 0.2e1 * t7 * t19 + 0.2e1 * t163 * t20 + 0.2e1 * t6 * t21 + 0.2e1 * t17 * t15 + 0.2e1 * t13 * t18 + 0.2e1 * t69 * t114 + (t120 * t69 + t24 * t36 + t148) * t158 + 0.2e1 * m(6) * (t13 * t4 + t163 * t3 + t148) + 0.2e1 * m(7) * (t1 * t6 + t17 * t5 + t2 * t7) - (mrSges(5,3) * t154 + t134 * t91 + t135 * t89) * t58 + (mrSges(5,1) * t112 + t24 * t157 - (-Ifges(6,6) * t89 + t156) * t58 + ((2 * Ifges(5,2)) + t166) * t57 + t98) * t61 + (mrSges(5,2) * t112 + mrSges(5,3) * t155 - 0.2e1 * Ifges(5,1) * t58 + (t10 + t11) * t91 + (t8 - t9) * t89 + (t156 + t138 * t91 + (-Ifges(6,6) + Ifges(7,6)) * t89) * t57 + ((t135 - t144) * t91 + (-t138 * t61 - t134) * t89) * qJD(5)) * t62 + 0.2e1 * (t110 * t115 + (-t90 ^ 2 + t92 ^ 2) * Ifges(4,4) + (Ifges(4,1) - Ifges(4,2)) * t90 * t92) * qJD(3); (t15 + t16) * t61 + (t37 + t38) * t57 - t95 * t58 + m(7) * (-t6 * t139 - t7 * t141 + t17 * t57 + t5 * t61) + m(5) * (-t36 * t58 + t101) + m(6) * (t13 * t141 - t139 * t163 + t101) + (m(7) * (t125 * t7 - t126 * t6) + m(5) * t24 + (m(7) * t1 + t160) * t91 + (t159 + t164) * t89) * t62; (-t146 + t147) * t158 + 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (-t128 * t146 + t147); -t24 * mrSges(5,2) + m(7) * (-t17 * t56 + t5 * t59) - t56 * t37 - Ifges(5,6) * t57 - Ifges(5,5) * t58 + t59 * t15 + t17 * t63 + t35 * t64 + t5 * t70 + t79 * t16 + (m(6) * t79 + t130) * t23 + (Ifges(4,5) * t92 - Ifges(4,6) * t90 + (-mrSges(4,1) * t92 + mrSges(4,2) * t90) * t78) * qJD(3) + (m(5) * (-t23 * t87 + t24 * t86) + (-t86 * t57 + t87 * t58) * mrSges(5,3)) * pkin(3) + (-t4 * mrSges(6,3) + t2 * mrSges(7,2) + t10 / 0.2e1 + t11 / 0.2e1 + (t65 / 0.2e1 - t66 / 0.2e1) * t62 - t117 * t58 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t57 + (t25 / 0.2e1 - t26 / 0.2e1 - t163 * mrSges(6,3) - t6 * mrSges(7,2) - t144 / 0.2e1 - t116 * t62) * qJD(5) + (m(7) * (-t6 * qJD(5) + t2) + t159) * t77) * t89 + (t3 * mrSges(6,3) + t1 * mrSges(7,2) - t8 / 0.2e1 + t9 / 0.2e1 + (t67 / 0.2e1 + t68 / 0.2e1) * t62 - t116 * t58 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t57 + (t27 / 0.2e1 + t28 / 0.2e1 + t7 * mrSges(7,2) - t13 * mrSges(6,3) + t117 * t62) * qJD(5) + (m(7) * (t7 * qJD(5) + t1) + t160) * t77) * t91 + t162 * t61 / 0.2e1; t53 + (t63 + t64) * t61 - t110 * qJD(3) + (t70 + t130) * t57 + m(7) * (-t56 * t61 + t57 * t59 - t131) + m(6) * (t57 * t79 - t131) + m(5) * (-t57 * t87 - t58 * t86) * pkin(3) - (mrSges(7,2) + mrSges(6,3)) * t167; 0.2e1 * t59 * t63 + 0.2e1 * t64 * t79 + (-t65 + t66) * t91 + (t67 + t68) * t89 + 0.2e1 * (-m(7) * t59 - t70) * t56 + ((t74 + t75) * t91 + (t72 - t73) * t89) * qJD(5); m(5) * t120 + t137 * t91 + t136 * t89 + t95 * qJD(5) + m(7) * (t1 * t89 - t2 * t91 + (t6 * t91 + t7 * t89) * qJD(5)) + m(6) * (t3 * t89 + t4 * t91 + (-t13 * t89 + t163 * t91) * qJD(5)) + t114; 0; 0; 0; Ifges(6,6) * t141 - pkin(5) * t19 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) + qJD(6) * t42 + qJ(6) * t21 + t1 * mrSges(7,3) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t4 * mrSges(6,1) + (-Ifges(6,6) * t91 - t138 * t89) * t62 * qJD(5) + t98; -(-m(7) * t102 - t108 - t109) * t58 + t93 * t62; -t161 * mrSges(7,2) - Ifges(6,6) * t126 + t93 * t77 + t162; m(7) * t56 + ((-mrSges(6,2) + mrSges(7,3)) * t91 + (-mrSges(6,1) - mrSges(7,1)) * t89) * qJD(5); 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); t19 + t164; t97 * m(7); (m(7) * t77 + mrSges(7,2)) * t125; m(7) * t126; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;

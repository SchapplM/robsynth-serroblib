% Calculate time derivative of joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:02
% EndTime: 2019-12-31 21:13:07
% DurationCPUTime: 1.72s
% Computational Cost: add. (3603->240), mult. (8075->373), div. (0->0), fcn. (7600->8), ass. (0->121)
t101 = cos(qJ(5));
t98 = sin(qJ(5));
t82 = -mrSges(6,1) * t101 + mrSges(6,2) * t98;
t176 = -mrSges(5,1) + t82;
t131 = qJD(5) * t101;
t137 = cos(pkin(9));
t173 = qJD(2) + qJD(3);
t100 = sin(qJ(2));
t102 = cos(qJ(3));
t103 = cos(qJ(2));
t133 = t102 * t103;
t99 = sin(qJ(3));
t75 = -t99 * t100 + t133;
t57 = t173 * t75;
t76 = t102 * t100 + t103 * t99;
t58 = t173 * t76;
t97 = sin(pkin(9));
t40 = t137 * t57 - t97 * t58;
t54 = t137 * t76 + t97 * t75;
t109 = t54 * t131 + t40 * t98;
t163 = -pkin(7) - pkin(6);
t86 = t163 * t100;
t77 = t99 * t86;
t87 = t163 * t103;
t60 = -t102 * t87 + t77;
t145 = t101 ^ 2 + t98 ^ 2;
t143 = pkin(2) * qJD(3);
t149 = t97 * t99;
t65 = (t102 * t137 - t149) * t143;
t122 = t145 * t65;
t53 = -t137 * t75 + t76 * t97;
t92 = -pkin(2) * t103 - pkin(1);
t61 = -t75 * pkin(3) + t92;
t25 = t53 * pkin(4) - t54 * pkin(8) + t61;
t59 = t102 * t86 + t87 * t99;
t107 = -qJ(4) * t76 + t59;
t50 = qJ(4) * t75 + t60;
t27 = t107 * t97 + t137 * t50;
t13 = t101 * t25 - t27 * t98;
t135 = qJD(5) * t98;
t14 = t101 * t27 + t25 * t98;
t175 = -t13 * t131 - t14 * t135;
t174 = t101 * t14 - t13 * t98;
t129 = t54 * t135;
t140 = t101 * t40;
t110 = t129 - t140;
t39 = t137 * t58 + t57 * t97;
t15 = mrSges(6,1) * t39 + mrSges(6,3) * t110;
t16 = -mrSges(6,2) * t39 - mrSges(6,3) * t109;
t148 = t98 * mrSges(6,3);
t32 = -mrSges(6,2) * t53 - t148 * t54;
t139 = t101 * t54;
t33 = mrSges(6,1) * t53 - mrSges(6,3) * t139;
t172 = t101 * t16 - t33 * t131 - t32 * t135 - t98 * t15;
t171 = 2 * m(5);
t170 = 2 * m(6);
t111 = -t57 * qJ(4) - t76 * qJD(4);
t132 = qJD(3) * t102;
t136 = qJD(3) * t99;
t42 = t76 * qJD(2) * t163 + t86 * t132 + t136 * t87;
t23 = -qJ(4) * t58 + qJD(4) * t75 + t42;
t105 = (t133 * t163 - t77) * qJD(2);
t43 = -t60 * qJD(3) + t105;
t11 = t23 * t97 - t137 * (t111 + t43);
t169 = 0.2e1 * t11;
t26 = -t107 * t137 + t50 * t97;
t168 = 0.2e1 * t26;
t51 = qJD(2) * t100 * pkin(2) + pkin(3) * t58;
t167 = 0.2e1 * t51;
t166 = 0.2e1 * t92;
t165 = m(5) * pkin(3);
t162 = pkin(3) * t97;
t12 = t137 * t23 + (t132 * t87 - t136 * t86 + t105 + t111) * t97;
t17 = pkin(4) * t39 - pkin(8) * t40 + t51;
t3 = -qJD(5) * t14 + t101 * t17 - t12 * t98;
t161 = t3 * t98;
t160 = Ifges(6,4) * t98;
t159 = Ifges(6,6) * t98;
t2 = qJD(5) * t13 + t101 * t12 + t17 * t98;
t158 = t101 * t2;
t157 = t11 * t26;
t117 = t137 * t99;
t64 = (t102 * t97 + t117) * t143;
t155 = t26 * t64;
t154 = t39 * mrSges(5,3);
t91 = pkin(2) * t102 + pkin(3);
t66 = -pkin(2) * t149 + t137 * t91;
t62 = -pkin(4) - t66;
t116 = mrSges(6,1) * t98 + mrSges(6,2) * t101;
t79 = t116 * qJD(5);
t152 = t62 * t79;
t121 = t137 * pkin(3);
t90 = -t121 - pkin(4);
t150 = t90 * t79;
t146 = Ifges(6,5) * t140 + Ifges(6,3) * t39;
t67 = pkin(2) * t117 + t97 * t91;
t144 = Ifges(6,4) * t101;
t130 = 0.2e1 * t103;
t119 = t39 * mrSges(5,1) + t40 * mrSges(5,2);
t118 = -(2 * Ifges(5,4)) - t159;
t115 = Ifges(6,1) * t101 - t160;
t114 = -Ifges(6,2) * t98 + t144;
t113 = Ifges(6,5) * t98 + Ifges(6,6) * t101;
t112 = t101 * t32 - t98 * t33;
t80 = t114 * qJD(5);
t81 = t115 * qJD(5);
t83 = Ifges(6,2) * t101 + t160;
t84 = Ifges(6,1) * t98 + t144;
t108 = t101 * t80 + t84 * t131 - t135 * t83 + t98 * t81;
t106 = -t161 + (-t101 * t13 - t14 * t98) * qJD(5);
t20 = t53 * Ifges(6,6) + t114 * t54;
t21 = t53 * Ifges(6,5) + t115 * t54;
t6 = -Ifges(6,4) * t110 - Ifges(6,2) * t109 + t39 * Ifges(6,6);
t7 = -Ifges(6,1) * t110 - Ifges(6,4) * t109 + t39 * Ifges(6,5);
t93 = Ifges(6,5) * t131;
t104 = -t42 * mrSges(4,2) - t12 * mrSges(5,2) + mrSges(6,3) * t158 + t21 * t131 / 0.2e1 + t26 * t79 + t84 * t140 / 0.2e1 + Ifges(5,5) * t40 + t43 * mrSges(4,1) + t81 * t139 / 0.2e1 + t53 * (-Ifges(6,6) * t135 + t93) / 0.2e1 + t101 * t6 / 0.2e1 - Ifges(4,6) * t58 + Ifges(4,5) * t57 + (t7 / 0.2e1 - t54 * t80 / 0.2e1) * t98 + (t113 / 0.2e1 - Ifges(5,6)) * t39 - t109 * t83 / 0.2e1 - (t54 * t84 + t20) * t135 / 0.2e1 + t176 * t11;
t89 = pkin(8) + t162;
t63 = pkin(8) + t67;
t31 = t116 * t54;
t9 = mrSges(6,1) * t109 - mrSges(6,2) * t110;
t1 = [-0.2e1 * t27 * t154 + 0.2e1 * t13 * t15 + 0.2e1 * t14 * t16 + t9 * t168 + t31 * t169 + 0.2e1 * t2 * t32 + 0.2e1 * t3 * t33 + 0.2e1 * t61 * t119 - 0.2e1 * t75 * Ifges(4,2) * t58 + 0.2e1 * t76 * t57 * Ifges(4,1) + (t58 * mrSges(4,1) + t57 * mrSges(4,2)) * t166 + (mrSges(5,3) * t168 + t101 * t21 - t98 * t20) * t40 + (t12 * t27 + t51 * t61 + t157) * t171 + 0.2e1 * m(4) * (t42 * t60 + t43 * t59) + (t13 * t3 + t14 * t2 + t157) * t170 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t103) * t130 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t75 + mrSges(4,2) * t76) + m(4) * pkin(2) * t166 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t100 + (-Ifges(3,2) + Ifges(3,1)) * t130) * t100) * qJD(2) + (mrSges(5,1) * t167 - 0.2e1 * mrSges(5,3) * t12 + t118 * t40 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t39 + t146) * t53 + (mrSges(5,2) * t167 + mrSges(5,3) * t169 + 0.2e1 * Ifges(5,1) * t40 + t101 * t7 - t98 * t6 + (Ifges(6,5) * t101 + t118) * t39 + (-t101 * t20 - t113 * t53 - t98 * t21) * qJD(5)) * t54 + 0.2e1 * (t75 * t57 - t76 * t58) * Ifges(4,4) + 0.2e1 * (t42 * t75 - t43 * t76 - t59 * t57 - t60 * t58) * mrSges(4,3); (Ifges(3,5) * t103 - Ifges(3,6) * t100 + (-mrSges(3,1) * t103 + mrSges(3,2) * t100) * pkin(6)) * qJD(2) + m(5) * (-t11 * t66 + t12 * t67 + t155) + m(6) * (t11 * t62 + t155) + (m(4) * (t102 * t43 + t42 * t99 + (t102 * t60 - t59 * t99) * qJD(3)) + (-t102 * t57 - t99 * t58 + (t102 * t75 + t76 * t99) * qJD(3)) * mrSges(4,3)) * pkin(2) + (m(6) * (t158 - t161 + t175) + t172) * t63 + t104 + (-t67 * t39 - t66 * t40 + t64 * t54) * mrSges(5,3) + t106 * mrSges(6,3) + t62 * t9 + t64 * t31 + (m(5) * t27 + m(6) * t174 - t53 * mrSges(5,3) + t112) * t65; 0.2e1 * t152 - 0.2e1 * t65 * mrSges(5,2) + (-t64 * t66 + t65 * t67) * t171 + (t122 * t63 + t62 * t64) * t170 + t108 + 0.2e1 * t176 * t64 + 0.2e1 * (-mrSges(4,1) * t99 - mrSges(4,2) * t102) * t143 + 0.2e1 * mrSges(6,3) * t122; -t154 * t162 - t40 * mrSges(5,3) * t121 + t104 + (-t11 * t137 + t12 * t97) * t165 - t3 * t148 + (m(6) * t11 + t9) * t90 + t175 * mrSges(6,3) + (m(6) * (t106 + t158) + t172) * t89; m(6) * t89 * t122 + t108 + t150 + t152 + (m(6) * t90 - t137 * t165 + t176) * t64 + (-mrSges(4,1) * t136 - mrSges(4,2) * t132) * pkin(2) + (mrSges(6,3) * t145 + t165 * t97 - mrSges(5,2)) * t65; t108 + 0.2e1 * t150; t101 * t15 + t98 * t16 + t112 * qJD(5) + m(6) * (t174 * qJD(5) + t101 * t3 + t2 * t98) + m(5) * t51 + t119; 0; 0; 0; mrSges(6,1) * t3 - mrSges(6,2) * t2 - Ifges(6,5) * t129 - Ifges(6,6) * t109 + t146; t93 - t116 * t65 + (t63 * t82 - t159) * qJD(5); t93 + (t82 * t89 - t159) * qJD(5); -t79; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

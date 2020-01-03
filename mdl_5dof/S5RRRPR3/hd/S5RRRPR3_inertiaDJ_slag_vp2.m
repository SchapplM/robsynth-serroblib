% Calculate time derivative of joint inertia matrix for
% S5RRRPR3
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:51
% EndTime: 2020-01-03 12:08:55
% DurationCPUTime: 1.14s
% Computational Cost: add. (2336->213), mult. (5086->291), div. (0->0), fcn. (4445->8), ass. (0->112)
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t157 = t108 ^ 2 + t111 ^ 2;
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t105 = sin(pkin(9));
t106 = cos(pkin(9));
t83 = t105 * t111 + t106 * t108;
t145 = pkin(8) * t83;
t109 = sin(qJ(2));
t95 = pkin(1) * t109 + pkin(7);
t129 = -qJ(4) - t95;
t79 = t129 * t108;
t102 = t111 * qJ(4);
t131 = t111 * t95;
t80 = t102 + t131;
t40 = -t105 * t80 + t106 * t79;
t30 = t40 - t145;
t41 = t105 * t79 + t106 * t80;
t82 = -t105 * t108 + t106 * t111;
t78 = t82 * pkin(8);
t31 = t78 + t41;
t12 = -t107 * t31 + t110 * t30;
t75 = t82 * qJD(3);
t146 = pkin(8) * t75;
t101 = t111 * qJD(4);
t119 = qJD(3) * t129;
t112 = cos(qJ(2));
t134 = pkin(1) * qJD(2);
t121 = t112 * t134;
t49 = t108 * t119 + t111 * t121 + t101;
t50 = (-qJD(4) - t121) * t108 + t111 * t119;
t25 = -t105 * t49 + t106 * t50;
t17 = t25 - t146;
t26 = t105 * t50 + t106 * t49;
t74 = t83 * qJD(3);
t69 = t74 * pkin(8);
t18 = -t69 + t26;
t2 = t12 * qJD(5) + t107 * t17 + t110 * t18;
t13 = t107 * t30 + t110 * t31;
t3 = -t13 * qJD(5) - t107 * t18 + t110 * t17;
t156 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t138 = -qJ(4) - pkin(7);
t90 = t138 * t108;
t92 = pkin(7) * t111 + t102;
t52 = -t105 * t92 + t106 * t90;
t38 = t52 - t145;
t53 = t105 * t90 + t106 * t92;
t39 = t78 + t53;
t15 = -t107 * t39 + t110 * t38;
t120 = qJD(3) * t138;
t72 = t108 * t120 + t101;
t73 = -t108 * qJD(4) + t111 * t120;
t36 = -t105 * t72 + t106 * t73;
t27 = t36 - t146;
t37 = t105 * t73 + t106 * t72;
t28 = -t69 + t37;
t5 = t15 * qJD(5) + t107 * t27 + t110 * t28;
t16 = t107 * t38 + t110 * t39;
t6 = -t16 * qJD(5) - t107 * t28 + t110 * t27;
t155 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t154 = t157 * t112;
t91 = -mrSges(4,1) * t111 + mrSges(4,2) * t108;
t153 = (t157 * mrSges(4,3) - mrSges(3,2)) * t112 + (t91 - mrSges(3,1)) * t109;
t152 = 2 * m(5);
t151 = 2 * m(6);
t43 = -t107 * t83 + t110 * t82;
t22 = t43 * qJD(5) - t107 * t74 + t110 * t75;
t44 = t107 * t82 + t110 * t83;
t23 = -t44 * qJD(5) - t107 * t75 - t110 * t74;
t9 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t150 = 0.2e1 * t9;
t24 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t149 = 0.2e1 * t24;
t42 = t74 * mrSges(5,1) + t75 * mrSges(5,2);
t148 = 0.2e1 * t42;
t48 = -mrSges(5,1) * t82 + mrSges(5,2) * t83;
t147 = pkin(3) * t48;
t142 = t82 * pkin(4);
t141 = pkin(1) * t112;
t140 = pkin(3) * t105;
t100 = t109 * t134;
t126 = qJD(3) * t108;
t99 = pkin(3) * t126;
t88 = t100 + t99;
t139 = t88 * t48;
t137 = Ifges(6,5) * t22 + Ifges(6,6) * t23;
t136 = Ifges(4,4) * t108;
t125 = qJD(3) * t111;
t124 = 2 * mrSges(5,3);
t123 = 2 * mrSges(6,3);
t122 = t106 * t75 * mrSges(5,3);
t97 = -t111 * pkin(3) - pkin(2);
t55 = pkin(4) * t74 + t99;
t94 = pkin(3) * t106 + pkin(4);
t70 = -t107 * t140 + t110 * t94;
t61 = t70 * qJD(5);
t71 = t107 * t94 + t110 * t140;
t62 = t71 * qJD(5);
t117 = -t62 * mrSges(6,1) - t61 * mrSges(6,2);
t116 = mrSges(4,1) * t108 + mrSges(4,2) * t111;
t89 = t97 - t141;
t115 = t42 + t9;
t114 = 0.2e1 * t22 * t44 * Ifges(6,1) + 0.2e1 * t23 * Ifges(6,2) * t43 - 0.2e1 * t82 * Ifges(5,2) * t74 + 0.2e1 * t83 * t75 * Ifges(5,1) + (Ifges(4,1) * t111 - t136) * t126 + (0.2e1 * Ifges(4,4) * t111 + (Ifges(4,1) - Ifges(4,2)) * t108) * t125 + 0.2e1 * (t22 * t43 + t23 * t44) * Ifges(6,4) + 0.2e1 * (-t83 * t74 + t82 * t75) * Ifges(5,4);
t113 = Ifges(4,5) * t125 + Ifges(5,5) * t75 + t137 + (-mrSges(5,3) * t140 - Ifges(5,6)) * t74 + (-t70 * t22 + t23 * t71 + t43 * t61 + t62 * t44) * mrSges(6,3);
t96 = -pkin(2) - t141;
t93 = Ifges(4,2) * t111 + t136;
t87 = t116 * qJD(3);
t57 = t97 - t142;
t54 = t89 - t142;
t51 = t100 + t55;
t1 = [0.2e1 * t139 + t89 * t148 + 0.2e1 * t96 * t87 + t51 * t149 + t54 * t150 - t93 * t126 + t114 + 0.2e1 * (m(4) * (t109 * t96 + t154 * t95) + t153) * t134 + (-t12 * t22 + t13 * t23 + t2 * t43 - t3 * t44) * t123 + (-t25 * t83 + t26 * t82 - t40 * t75 - t41 * t74) * t124 + (t25 * t40 + t26 * t41 + t88 * t89) * t152 + (t12 * t3 + t13 * t2 + t51 * t54) * t151; t139 + t114 + m(5) * (t25 * t52 + t26 * t53 + t36 * t40 + t37 * t41 + t88 * t97 + t89 * t99) + m(6) * (t12 * t6 + t13 * t5 + t15 * t3 + t16 * t2 + t51 * t57 + t54 * t55) + ((-t25 - t36) * t83 + (t26 + t37) * t82 + (-t40 - t52) * t75 - (t41 + t53) * t74) * mrSges(5,3) + ((-t3 - t6) * t44 + (t2 + t5) * t43 + (t13 + t16) * t23 + (-t12 - t15) * t22) * mrSges(6,3) + (t57 + t54) * t9 + (-pkin(2) + t96) * t87 + (t89 + t97) * t42 + (t51 + t55) * t24 + (m(4) * (-pkin(2) * t109 + pkin(7) * t154) + t153) * t134 + (-t93 + t147) * t126; -0.2e1 * pkin(2) * t87 + t97 * t148 + t57 * t150 + t55 * t149 + t114 + (-t15 * t22 + t16 * t23 + t5 * t43 - t6 * t44) * t123 + (-t36 * t83 + t37 * t82 - t52 * t75 - t53 * t74) * t124 + (-t93 + 0.2e1 * t147) * t126 + (t36 * t52 + t37 * t53 + t97 * t99) * t152 + (t15 * t6 + t16 * t5 + t55 * t57) * t151; m(6) * (-t12 * t62 + t13 * t61 + t2 * t71 + t3 * t70) + (-t122 + m(5) * (t105 * t26 + t106 * t25)) * pkin(3) + t25 * mrSges(5,1) - t26 * mrSges(5,2) + t113 + (-mrSges(4,1) * t131 + (mrSges(4,2) * t95 - Ifges(4,6)) * t108) * qJD(3) - t116 * t121 + t156; (-Ifges(4,6) * t108 + t91 * pkin(7)) * qJD(3) + m(6) * (-t15 * t62 + t16 * t61 + t5 * t71 + t6 * t70) + (-t122 + m(5) * (t105 * t37 + t106 * t36)) * pkin(3) + t36 * mrSges(5,1) - t37 * mrSges(5,2) + t113 + t155; 0.2e1 * m(6) * (t61 * t71 - t62 * t70) + 0.2e1 * t117; m(5) * t88 + m(6) * t51 + t115; m(5) * t99 + m(6) * t55 + t115; 0; 0; t137 + t156; t137 + t155; t117; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

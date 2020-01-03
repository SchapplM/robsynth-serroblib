% Calculate time derivative of joint inertia matrix for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:04
% EndTime: 2019-12-31 20:20:09
% DurationCPUTime: 2.03s
% Computational Cost: add. (3336->297), mult. (7404->460), div. (0->0), fcn. (6918->8), ass. (0->132)
t107 = cos(qJ(4));
t101 = sin(pkin(9));
t102 = cos(pkin(9));
t105 = sin(qJ(2));
t108 = cos(qJ(2));
t83 = t101 * t105 - t102 * t108;
t79 = t83 * qJD(2);
t135 = t107 * t79;
t84 = t101 * t108 + t102 * t105;
t78 = t84 * qJD(2);
t157 = -Ifges(5,5) * t135 + Ifges(5,3) * t78;
t103 = sin(qJ(5));
t104 = sin(qJ(4));
t106 = cos(qJ(5));
t113 = t103 * t104 - t106 * t107;
t51 = t113 * t84;
t97 = pkin(2) * t101 + pkin(7);
t156 = -m(5) * t97 - mrSges(5,3);
t117 = mrSges(5,1) * t104 + mrSges(5,2) * t107;
t88 = t117 * qJD(4);
t155 = qJD(4) + qJD(5);
t154 = 2 * m(6);
t153 = -2 * mrSges(4,3);
t140 = -qJ(3) - pkin(6);
t92 = t140 * t105;
t93 = t140 * t108;
t68 = -t101 * t93 - t102 * t92;
t151 = 0.2e1 * t68;
t150 = m(6) * pkin(4);
t148 = -t84 / 0.2e1;
t138 = Ifges(5,4) * t104;
t94 = Ifges(5,2) * t107 + t138;
t147 = -t94 / 0.2e1;
t145 = pkin(8) + t97;
t119 = qJD(2) * t140;
t109 = -t105 * qJD(3) + t108 * t119;
t77 = qJD(3) * t108 + t105 * t119;
t46 = t101 * t77 - t102 * t109;
t144 = t46 * t68;
t143 = t78 * Ifges(5,5);
t142 = t78 * Ifges(5,6);
t141 = t83 * Ifges(5,6);
t63 = t155 * t113;
t87 = t103 * t107 + t104 * t106;
t64 = t155 * t87;
t139 = -Ifges(6,5) * t63 - Ifges(6,6) * t64;
t99 = -pkin(2) * t108 - pkin(1);
t56 = t83 * pkin(3) - t84 * pkin(7) + t99;
t69 = t101 * t92 - t102 * t93;
t62 = t107 * t69;
t32 = t104 * t56 + t62;
t137 = Ifges(5,4) * t107;
t136 = t104 * t84;
t134 = t107 * t84;
t132 = qJD(4) * t104;
t131 = qJD(4) * t107;
t130 = qJD(5) * t103;
t129 = qJD(5) * t106;
t17 = t113 * t79 - t64 * t84;
t18 = t155 * t51 + t87 * t79;
t128 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t78;
t126 = pkin(2) * qJD(2) * t105;
t125 = pkin(4) * t132;
t124 = t84 * t132;
t98 = -pkin(2) * t102 - pkin(3);
t123 = t78 * mrSges(4,1) - t79 * mrSges(4,2);
t33 = t64 * mrSges(6,1) - t63 * mrSges(6,2);
t122 = qJD(4) * t145;
t121 = -Ifges(5,6) * t104 - (2 * Ifges(4,4));
t47 = t101 * t109 + t102 * t77;
t48 = pkin(3) * t78 + pkin(7) * t79 + t126;
t120 = -t104 * t47 + t107 * t48;
t31 = -t104 * t69 + t107 * t56;
t118 = 0.2e1 * t126;
t116 = Ifges(5,1) * t107 - t138;
t115 = -Ifges(5,2) * t104 + t137;
t19 = pkin(4) * t83 - pkin(8) * t134 + t31;
t24 = -pkin(8) * t136 + t32;
t9 = -t103 * t24 + t106 * t19;
t10 = t103 * t19 + t106 * t24;
t80 = t145 * t104;
t81 = t145 * t107;
t54 = -t103 * t81 - t106 * t80;
t55 = -t103 * t80 + t106 * t81;
t75 = t104 * t122;
t76 = t107 * t122;
t29 = t54 * qJD(5) - t103 * t76 - t106 * t75;
t30 = -t55 * qJD(5) + t103 * t75 - t106 * t76;
t114 = t30 * mrSges(6,1) - t29 * mrSges(6,2) + t139;
t7 = pkin(8) * t135 + pkin(4) * t78 + (-t62 + (pkin(8) * t84 - t56) * t104) * qJD(4) + t120;
t11 = t104 * t48 + t107 * t47 + t56 * t131 - t69 * t132;
t111 = -t104 * t79 + t84 * t131;
t8 = -t111 * pkin(8) + t11;
t2 = t9 * qJD(5) + t103 * t7 + t106 * t8;
t3 = -t10 * qJD(5) - t103 * t8 + t106 * t7;
t112 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t128;
t110 = t124 + t135;
t100 = Ifges(5,5) * t131;
t95 = Ifges(5,1) * t104 + t137;
t91 = -pkin(4) * t107 + t98;
t90 = t116 * qJD(4);
t89 = t115 * qJD(4);
t82 = (-mrSges(6,1) * t103 - mrSges(6,2) * t106) * qJD(5) * pkin(4);
t67 = Ifges(6,1) * t87 - Ifges(6,4) * t113;
t66 = Ifges(6,4) * t87 - Ifges(6,2) * t113;
t65 = mrSges(6,1) * t113 + mrSges(6,2) * t87;
t58 = mrSges(5,1) * t83 - mrSges(5,3) * t134;
t57 = -mrSges(5,2) * t83 - mrSges(5,3) * t136;
t50 = t87 * t84;
t45 = pkin(4) * t136 + t68;
t41 = Ifges(5,5) * t83 + t116 * t84;
t40 = t115 * t84 + t141;
t39 = mrSges(6,1) * t83 + mrSges(6,3) * t51;
t38 = -mrSges(6,2) * t83 - mrSges(6,3) * t50;
t37 = -mrSges(5,2) * t78 - t111 * mrSges(5,3);
t36 = mrSges(5,1) * t78 + t110 * mrSges(5,3);
t35 = -Ifges(6,1) * t63 - Ifges(6,4) * t64;
t34 = -Ifges(6,4) * t63 - Ifges(6,2) * t64;
t28 = t111 * mrSges(5,1) - t110 * mrSges(5,2);
t27 = t111 * pkin(4) + t46;
t25 = mrSges(6,1) * t50 - mrSges(6,2) * t51;
t23 = -Ifges(6,1) * t51 - Ifges(6,4) * t50 + Ifges(6,5) * t83;
t22 = -Ifges(6,4) * t51 - Ifges(6,2) * t50 + Ifges(6,6) * t83;
t21 = -t110 * Ifges(5,1) - t111 * Ifges(5,4) + t143;
t20 = -t110 * Ifges(5,4) - t111 * Ifges(5,2) + t142;
t14 = -mrSges(6,2) * t78 + mrSges(6,3) * t18;
t13 = mrSges(6,1) * t78 - mrSges(6,3) * t17;
t12 = -t32 * qJD(4) + t120;
t6 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t5 = Ifges(6,1) * t17 + Ifges(6,4) * t18 + t78 * Ifges(6,5);
t4 = Ifges(6,4) * t17 + Ifges(6,2) * t18 + t78 * Ifges(6,6);
t1 = [(mrSges(4,2) * t118 - t104 * t20 + t107 * t21 - 0.2e1 * Ifges(4,1) * t79 + (-t107 * t40 - t104 * t41 + t83 * (-Ifges(5,5) * t104 - Ifges(5,6) * t107)) * qJD(4) + 0.2e1 * (mrSges(4,3) + t117) * t46) * t84 + 0.2e1 * t99 * t123 + t28 * t151 + (t10 * t2 + t27 * t45 + t3 * t9) * t154 + 0.2e1 * t11 * t57 + 0.2e1 * t12 * t58 + (mrSges(4,1) * t118 - t121 * t79 + t47 * t153 + t128 + t157) * t83 - t50 * t4 - t51 * t5 + 0.2e1 * t31 * t36 + 0.2e1 * t32 * t37 + 0.2e1 * t2 * t38 + 0.2e1 * t3 * t39 + 0.2e1 * t45 * t6 + 0.2e1 * t27 * t25 + t18 * t22 + t17 * t23 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + 0.2e1 * m(4) * (t99 * t126 + t47 * t69 + t144) + 0.2e1 * m(5) * (t11 * t32 + t12 * t31 + t144) - (mrSges(4,3) * t151 - t104 * t40 + t107 * t41) * t79 + (-Ifges(6,5) * t51 - Ifges(6,6) * t50 + (Ifges(5,5) * t107 + t121) * t84 + t69 * t153 + ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t83) * t78 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t105 + mrSges(3,2) * t108) + (-Ifges(3,2) + Ifges(3,1)) * t105 * t108 + (-t105 ^ 2 + t108 ^ 2) * Ifges(3,4)) * qJD(2); (t139 + t100) * t83 / 0.2e1 + t78 * (Ifges(6,5) * t87 - Ifges(6,6) * t113) / 0.2e1 - t113 * t4 / 0.2e1 + (-t10 * t64 - t113 * t2 - t3 * t87 + t9 * t63) * mrSges(6,3) + (m(4) * (t101 * t47 - t102 * t46) + (-t101 * t78 + t102 * t79) * mrSges(4,3)) * pkin(2) + (t11 * mrSges(5,3) + t84 * t90 / 0.2e1 - t79 * t95 / 0.2e1 + t142 / 0.2e1 - t46 * mrSges(5,1) + t20 / 0.2e1 + (t84 * t147 - t31 * mrSges(5,3) + t41 / 0.2e1) * qJD(4) + (-qJD(4) * t58 + m(5) * (-t31 * qJD(4) + t11) + t37) * t97) * t107 + (Ifges(3,5) * t108 - Ifges(3,6) * t105 + (-mrSges(3,1) * t108 + mrSges(3,2) * t105) * pkin(6)) * qJD(2) + m(6) * (t10 * t29 + t2 * t55 + t27 * t91 + t3 * t54 + t30 * t9) + t55 * t14 + (t89 * t148 - t79 * t147 - t97 * t36 + t21 / 0.2e1 + t143 / 0.2e1 + t46 * mrSges(5,2) + t156 * t12 + (t95 * t148 - t97 * t57 + pkin(4) * t25 - t40 / 0.2e1 - t141 / 0.2e1 + t45 * t150 + t156 * t32) * qJD(4)) * t104 - t63 * t23 / 0.2e1 - t64 * t22 / 0.2e1 + t27 * t65 + t18 * t66 / 0.2e1 + t17 * t67 / 0.2e1 - Ifges(4,6) * t78 - Ifges(4,5) * t79 + t87 * t5 / 0.2e1 + t68 * t88 + t91 * t6 - t46 * mrSges(4,1) - t47 * mrSges(4,2) - t50 * t34 / 0.2e1 - t51 * t35 / 0.2e1 + t54 * t13 + t29 * t38 + t30 * t39 + t45 * t33 + (m(5) * t46 + t28) * t98; -t63 * t67 + t87 * t35 - t64 * t66 - t113 * t34 + (t91 * t125 + t29 * t55 + t30 * t54) * t154 + 0.2e1 * t65 * t125 + 0.2e1 * t91 * t33 + 0.2e1 * t98 * t88 + t104 * t90 - t94 * t132 + (qJD(4) * t95 + t89) * t107 + 0.2e1 * (-t113 * t29 - t30 * t87 + t54 * t63 - t55 * t64) * mrSges(6,3); m(4) * t126 + t104 * t37 + t107 * t36 - t113 * t13 + t87 * t14 - t63 * t38 - t64 * t39 + (-t104 * t58 + t107 * t57) * qJD(4) + m(6) * (-t10 * t63 - t113 * t3 + t2 * t87 - t64 * t9) + m(5) * (t104 * t11 + t107 * t12 + (-t104 * t31 + t107 * t32) * qJD(4)) + t123; m(6) * (-t113 * t30 + t29 * t87 - t54 * t64 - t55 * t63); (t113 * t64 - t63 * t87) * t154; -Ifges(5,5) * t124 + t12 * mrSges(5,1) - t11 * mrSges(5,2) - t111 * Ifges(5,6) + (m(6) * (t10 * t129 + t103 * t2 + t106 * t3 - t9 * t130) + t38 * t129 + t103 * t14 - t39 * t130 + t106 * t13) * pkin(4) + t112 + t157; t100 + (-t107 * t97 * mrSges(5,1) + (mrSges(5,2) * t97 - Ifges(5,6)) * t104) * qJD(4) + (m(6) * (t103 * t29 + t106 * t30 + (-t103 * t54 + t106 * t55) * qJD(5)) + (-t103 * t64 + t106 * t63 + (t103 * t87 - t106 * t113) * qJD(5)) * mrSges(6,3)) * pkin(4) + t114; -t88 + (-t103 * t63 - t106 * t64 + (t103 * t113 + t106 * t87) * qJD(5)) * t150 - t33; 0.2e1 * t82; t112; t114; -t33; t82; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

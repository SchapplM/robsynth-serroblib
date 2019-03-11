% Calculate time derivative of joint inertia matrix for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:12
% EndTime: 2019-03-09 02:10:14
% DurationCPUTime: 1.65s
% Computational Cost: add. (1061->287), mult. (2267->400), div. (0->0), fcn. (1452->4), ass. (0->118)
t139 = Ifges(7,4) + Ifges(6,5);
t74 = sin(qJ(5));
t76 = cos(qJ(5));
t142 = -Ifges(7,6) * t74 - t139 * t76;
t141 = Ifges(7,2) + Ifges(6,3);
t114 = t74 ^ 2 + t76 ^ 2;
t77 = cos(qJ(4));
t109 = qJD(4) * t77;
t72 = qJ(2) - pkin(7);
t103 = t72 * t109;
t75 = sin(qJ(4));
t140 = qJD(2) * t75 + t103;
t86 = pkin(5) * t74 - qJ(6) * t76;
t26 = t86 * qJD(5) - qJD(6) * t74;
t107 = qJD(5) * t76;
t108 = qJD(5) * t74;
t33 = mrSges(7,1) * t108 - mrSges(7,3) * t107;
t34 = mrSges(6,1) * t108 + mrSges(6,2) * t107;
t138 = m(7) * t26 + t33 + t34;
t105 = qJD(6) * t76;
t87 = pkin(5) * t76 + qJ(6) * t74;
t137 = -t87 * qJD(5) + t105;
t125 = t72 * t75;
t102 = qJD(5) * t125;
t32 = qJD(3) + (pkin(4) * t77 + pkin(8) * t75) * qJD(4);
t73 = pkin(1) + qJ(3);
t43 = t75 * pkin(4) - pkin(8) * t77 + t73;
t4 = -(qJD(5) * t43 + t140) * t74 - (-t32 + t102) * t76;
t123 = t76 * t77;
t40 = t75 * mrSges(6,1) - mrSges(6,3) * t123;
t41 = -t75 * mrSges(7,1) + mrSges(7,2) * t123;
t117 = t40 - t41;
t124 = t74 * t77;
t39 = -t75 * mrSges(6,2) - mrSges(6,3) * t124;
t42 = -mrSges(7,2) * t124 + t75 * mrSges(7,3);
t118 = t39 + t42;
t136 = t117 * t74 - t118 * t76;
t135 = 0.2e1 * qJD(3);
t133 = Ifges(6,4) * t74;
t132 = Ifges(6,4) * t76;
t131 = Ifges(7,5) * t74;
t130 = Ifges(7,5) * t76;
t129 = Ifges(6,6) * t75;
t127 = Ifges(7,6) * t75;
t126 = t43 * t76;
t106 = qJD(5) * t77;
t110 = qJD(4) * t75;
t81 = t74 * t106 + t76 * t110;
t17 = mrSges(6,1) * t109 + t81 * mrSges(6,3);
t18 = -mrSges(7,1) * t109 - t81 * mrSges(7,2);
t121 = -t17 + t18;
t100 = t74 * t110;
t101 = t76 * t106;
t80 = t100 - t101;
t19 = -mrSges(6,2) * t109 + t80 * mrSges(6,3);
t20 = t80 * mrSges(7,2) + mrSges(7,3) * t109;
t120 = t19 + t20;
t88 = Ifges(7,3) * t74 + t130;
t22 = t88 * t77 + t127;
t89 = -Ifges(6,2) * t74 + t132;
t23 = t89 * t77 + t129;
t119 = t22 - t23;
t47 = -t76 * mrSges(6,1) + t74 * mrSges(6,2);
t116 = t47 - mrSges(5,1);
t16 = t76 * t125 + t74 * t43;
t115 = t114 * pkin(8) * t109;
t69 = t75 ^ 2;
t113 = qJD(2) * t69;
t111 = qJD(3) * t73;
t71 = t77 ^ 2;
t67 = t71 * qJD(2);
t48 = -Ifges(7,3) * t76 + t131;
t49 = Ifges(6,2) * t76 + t133;
t99 = t48 / 0.2e1 - t49 / 0.2e1;
t50 = Ifges(7,1) * t74 - t130;
t51 = Ifges(6,1) * t74 + t132;
t98 = -t50 / 0.2e1 - t51 / 0.2e1;
t97 = t43 * t107 + t140 * t76 + t74 * t32;
t96 = Ifges(6,6) * t100 + Ifges(7,6) * t101 + t141 * t109;
t95 = m(6) * pkin(4) - t116;
t93 = t74 * mrSges(6,1) + t76 * mrSges(6,2);
t46 = -t76 * mrSges(7,1) - t74 * mrSges(7,3);
t92 = t74 * mrSges(7,1) - t76 * mrSges(7,3);
t91 = Ifges(6,1) * t76 - t133;
t90 = Ifges(7,1) * t76 + t131;
t12 = qJ(6) * t75 + t16;
t13 = -t126 + (t72 * t74 - pkin(5)) * t75;
t85 = t12 * t76 + t13 * t74;
t15 = -t74 * t125 + t126;
t84 = t15 * t74 - t16 * t76;
t24 = Ifges(7,4) * t75 + t90 * t77;
t25 = Ifges(6,5) * t75 + t91 * t77;
t83 = -t139 * t75 - t24 - t25;
t82 = -t72 + t86;
t78 = m(7) * t105 + (-m(7) * t87 + t46 + t47) * qJD(5);
t66 = Ifges(7,4) * t107;
t65 = Ifges(6,5) * t107;
t63 = Ifges(7,6) * t108;
t57 = t72 * t67;
t44 = -pkin(4) - t87;
t38 = t91 * qJD(5);
t37 = t90 * qJD(5);
t36 = t89 * qJD(5);
t35 = t88 * qJD(5);
t30 = t93 * t77;
t29 = t92 * t77;
t21 = t82 * t77;
t11 = -t80 * mrSges(6,1) - t81 * mrSges(6,2);
t10 = -t80 * mrSges(7,1) + t81 * mrSges(7,3);
t9 = -t51 * t106 + (Ifges(6,5) * t77 - t91 * t75) * qJD(4);
t8 = -t50 * t106 + (Ifges(7,4) * t77 - t90 * t75) * qJD(4);
t7 = -t49 * t106 + (Ifges(6,6) * t77 - t89 * t75) * qJD(4);
t6 = -t48 * t106 + (Ifges(7,6) * t77 - t88 * t75) * qJD(4);
t5 = -t82 * t110 + (-qJD(2) - t137) * t77;
t3 = -t74 * t102 + t97;
t2 = -pkin(5) * t109 - t4;
t1 = qJ(6) * t109 + (-t72 * t108 + qJD(6)) * t75 + t97;
t14 = [mrSges(4,3) * t135 + 0.2e1 * t1 * t42 + 0.2e1 * t21 * t10 + 0.2e1 * t12 * t20 + 0.2e1 * t13 * t18 + 0.2e1 * t15 * t17 + 0.2e1 * t16 * t19 + 0.2e1 * t2 * t41 + 0.2e1 * t5 * t29 + 0.2e1 * t3 * t39 + 0.2e1 * t4 * t40 + 0.2e1 * (m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3) + (-t69 - t71) * mrSges(5,3)) * qJD(2) + 0.2e1 * m(5) * (t72 * t113 + t111 + t57) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t111) + 0.2e1 * m(7) * (t1 * t12 + t13 * t2 + t21 * t5) + 0.2e1 * m(6) * (t15 * t4 + t16 * t3 + t57) + (mrSges(5,1) * t135 + (-0.2e1 * t73 * mrSges(5,2) + 0.2e1 * Ifges(5,4) * t75 + 0.2e1 * t72 * t30 + (-t119 - t127) * t74 + t83 * t76) * qJD(4) + t96) * t75 + (mrSges(5,2) * t135 - 0.2e1 * qJD(2) * t30 - 0.2e1 * t72 * t11 + (t8 + t9) * t76 + (t6 - t7) * t74 + ((t119 - t129) * t76 + t83 * t74) * qJD(5) + (0.2e1 * t73 * mrSges(5,1) + (-Ifges(6,6) * t74 - 0.2e1 * Ifges(5,4) - t142) * t77 + (-0.2e1 * m(6) * t72 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + t141) * t75) * qJD(4)) * t77; t121 * t76 - t120 * t74 + (-mrSges(5,1) * t77 + mrSges(5,2) * t75) * qJD(4) + (-m(5) - m(4)) * qJD(3) + t136 * qJD(5) + m(7) * (-t85 * qJD(5) - t1 * t74 + t2 * t76) + m(6) * (t84 * qJD(5) - t3 * t74 - t4 * t76); 0; m(4) * qJD(2) + m(6) * t67 + m(5) * (t67 + t113) + (-m(7) * t5 - t11 - t10 + (-m(6) * t84 + m(7) * t85 - t136) * qJD(4)) * t77 + (t120 * t76 + t121 * t74 + (t29 + t30) * qJD(4) + (-t117 * t76 - t118 * t74) * qJD(5) + m(7) * (qJD(4) * t21 + t1 * t76 + t13 * t107 - t12 * t108 + t2 * t74) + m(6) * (-t15 * t107 - t16 * t108 + t3 * t76 - t4 * t74 - 0.2e1 * t103)) * t75; 0; 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (-0.1e1 + t114) * t75 * t109; t21 * t33 + t44 * t10 + t5 * t46 + t26 * t29 - pkin(4) * t11 + m(7) * (t21 * t26 + t44 * t5) + (-qJD(2) * mrSges(5,2) + t65 / 0.2e1 + t66 / 0.2e1 + t63 / 0.2e1 + (-t95 * t72 - Ifges(5,5)) * qJD(4)) * t75 + (-t4 * mrSges(6,3) + t2 * mrSges(7,2) + t8 / 0.2e1 + t9 / 0.2e1 - t99 * t110 + (-t12 * mrSges(7,2) - t16 * mrSges(6,3) + t22 / 0.2e1 - t23 / 0.2e1 - t129 / 0.2e1) * qJD(5) + (-t118 * qJD(5) + m(7) * (-qJD(5) * t12 + t2) + m(6) * (-qJD(5) * t16 - t4) + t121) * pkin(8)) * t74 + (t3 * mrSges(6,3) + t1 * mrSges(7,2) - t6 / 0.2e1 + t7 / 0.2e1 + t98 * t110 + (t13 * mrSges(7,2) - t15 * mrSges(6,3) + t25 / 0.2e1 + t24 / 0.2e1) * qJD(5) + (-t117 * qJD(5) + m(7) * (qJD(5) * t13 + t1) + m(6) * (-qJD(5) * t15 + t3) + t120) * pkin(8)) * t76 + (-t72 * t34 + (t37 / 0.2e1 + t38 / 0.2e1) * t76 + (t35 / 0.2e1 - t36 / 0.2e1) * t74 + (-t72 * mrSges(5,2) - Ifges(5,6) + (Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t76 + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t74) * qJD(4) + t95 * qJD(2) + (t98 * t74 + t99 * t76) * qJD(5)) * t77; 0; (t46 + t116) * t110 + m(7) * (t44 * t110 + t115) + m(6) * (-pkin(4) * t110 + t115) + ((-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t114) * qJD(4) - t138) * t77; -0.2e1 * pkin(4) * t34 + 0.2e1 * t33 * t44 + (-t35 + t36) * t76 + (t37 + t38) * t74 + 0.2e1 * (m(7) * t44 + t46) * t26 + ((t50 + t51) * t76 + (t48 - t49) * t74) * qJD(5); m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t12) - t3 * mrSges(6,2) - pkin(5) * t18 + t1 * mrSges(7,3) + qJD(6) * t42 + qJ(6) * t20 + t4 * mrSges(6,1) - t2 * mrSges(7,1) + (-Ifges(6,6) * t76 - t139 * t74) * t106 + t142 * t110 + t96; t138; (-m(7) * t86 - t92 - t93) * t109 + t78 * t75; t137 * mrSges(7,2) - Ifges(6,6) * t108 + t78 * pkin(8) + t63 + t65 + t66; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t2 + t18; -m(7) * t108; (t75 * t107 + t74 * t109) * m(7); (m(7) * pkin(8) + mrSges(7,2)) * t107; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;

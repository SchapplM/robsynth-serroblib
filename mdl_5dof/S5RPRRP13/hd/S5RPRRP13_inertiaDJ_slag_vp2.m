% Calculate time derivative of joint inertia matrix for
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:36
% DurationCPUTime: 1.47s
% Computational Cost: add. (844->242), mult. (1945->344), div. (0->0), fcn. (1243->4), ass. (0->109)
t130 = Ifges(6,4) + Ifges(5,5);
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t134 = -Ifges(6,6) * t64 - t130 * t66;
t124 = 2 * qJD(2);
t133 = Ifges(6,2) + Ifges(5,3);
t103 = t64 ^ 2 + t66 ^ 2;
t65 = sin(qJ(3));
t68 = -pkin(1) - pkin(6);
t115 = t65 * t68;
t67 = cos(qJ(3));
t43 = pkin(3) * t65 - pkin(7) * t67 + qJ(2);
t129 = t66 * t115 + t64 * t43;
t132 = qJD(4) * t129;
t131 = -2 * m(5);
t76 = pkin(4) * t66 + qJ(5) * t64;
t97 = qJD(5) * t66;
t128 = t76 * qJD(4) - t97;
t32 = qJD(2) + (pkin(3) * t67 + pkin(7) * t65) * qJD(3);
t127 = -t66 * t32 + t132;
t101 = qJD(3) * t67;
t100 = qJD(4) * t64;
t89 = t68 * t100;
t90 = t68 * t101;
t99 = qJD(4) * t66;
t92 = t64 * t32 + t43 * t99 + t66 * t90;
t1 = qJ(5) * t101 + (qJD(5) - t89) * t65 + t92;
t112 = t66 * t67;
t40 = mrSges(5,1) * t65 - mrSges(5,3) * t112;
t41 = -mrSges(6,1) * t65 + mrSges(6,2) * t112;
t106 = -t40 + t41;
t113 = t66 * t43;
t84 = t64 * t68 - pkin(4);
t14 = t84 * t65 - t113;
t102 = qJD(3) * t65;
t87 = t64 * t102;
t98 = qJD(4) * t67;
t88 = t66 * t98;
t70 = t87 - t88;
t17 = -mrSges(5,2) * t101 + t70 * mrSges(5,3);
t18 = t70 * mrSges(6,2) + mrSges(6,3) * t101;
t19 = -t64 * t115 + t113;
t3 = -t65 * t89 + t92;
t126 = t106 * qJD(4) + m(6) * (t14 * qJD(4) + t1) + m(5) * (-t19 * qJD(4) + t3) + t17 + t18;
t117 = t64 * t67;
t39 = -mrSges(5,2) * t65 - mrSges(5,3) * t117;
t42 = -mrSges(6,2) * t117 + mrSges(6,3) * t65;
t107 = t39 + t42;
t13 = qJ(5) * t65 + t129;
t71 = t66 * t102 + t64 * t98;
t15 = mrSges(5,1) * t101 + t71 * mrSges(5,3);
t16 = -mrSges(6,1) * t101 - t71 * mrSges(6,2);
t2 = t84 * t101 + t127;
t4 = -t64 * t90 - t127;
t125 = -t107 * qJD(4) + m(6) * (-t13 * qJD(4) + t2) + m(5) * (-t4 - t132) - t15 + t16;
t123 = Ifges(5,4) * t64;
t122 = Ifges(5,4) * t66;
t121 = Ifges(6,5) * t64;
t120 = Ifges(6,5) * t66;
t118 = Ifges(6,6) * t65;
t116 = t65 * Ifges(5,6);
t77 = Ifges(6,3) * t64 + t120;
t22 = t77 * t67 + t118;
t78 = -Ifges(5,2) * t64 + t122;
t23 = t78 * t67 + t116;
t108 = t22 - t23;
t47 = -t66 * mrSges(5,1) + t64 * mrSges(5,2);
t105 = t47 - mrSges(4,1);
t104 = t103 * pkin(7) * t101;
t48 = -Ifges(6,3) * t66 + t121;
t49 = Ifges(5,2) * t66 + t123;
t86 = t48 / 0.2e1 - t49 / 0.2e1;
t50 = Ifges(6,1) * t64 - t120;
t51 = Ifges(5,1) * t64 + t122;
t85 = -t50 / 0.2e1 - t51 / 0.2e1;
t83 = Ifges(5,6) * t87 + Ifges(6,6) * t88 + t133 * t101;
t82 = t64 * mrSges(5,1) + t66 * mrSges(5,2);
t46 = -t66 * mrSges(6,1) - t64 * mrSges(6,3);
t81 = t64 * mrSges(6,1) - t66 * mrSges(6,3);
t80 = Ifges(5,1) * t66 - t123;
t79 = Ifges(6,1) * t66 + t121;
t75 = pkin(4) * t64 - qJ(5) * t66;
t24 = t65 * Ifges(6,4) + t79 * t67;
t25 = t65 * Ifges(5,5) + t80 * t67;
t74 = -t130 * t65 - t24 - t25;
t72 = -t68 + t75;
t69 = m(6) * t97 + (-m(6) * t76 + t46 + t47) * qJD(4);
t61 = Ifges(6,4) * t99;
t60 = Ifges(5,5) * t99;
t58 = Ifges(6,6) * t100;
t44 = -pkin(3) - t76;
t38 = t80 * qJD(4);
t37 = t79 * qJD(4);
t36 = t78 * qJD(4);
t35 = t77 * qJD(4);
t34 = t82 * qJD(4);
t33 = t81 * qJD(4);
t30 = t82 * t67;
t29 = t81 * t67;
t26 = t75 * qJD(4) - qJD(5) * t64;
t21 = t72 * t67;
t11 = -t70 * mrSges(5,1) - t71 * mrSges(5,2);
t10 = -t70 * mrSges(6,1) + t71 * mrSges(6,3);
t9 = -t51 * t98 + (Ifges(5,5) * t67 - t80 * t65) * qJD(3);
t8 = -t50 * t98 + (Ifges(6,4) * t67 - t79 * t65) * qJD(3);
t7 = -t49 * t98 + (Ifges(5,6) * t67 - t78 * t65) * qJD(3);
t6 = -t48 * t98 + (Ifges(6,6) * t67 - t77 * t65) * qJD(3);
t5 = -t72 * t102 + t128 * t67;
t12 = [0.2e1 * t1 * t42 + 0.2e1 * t21 * t10 + 0.2e1 * t13 * t18 + 0.2e1 * t14 * t16 + 0.2e1 * t19 * t15 + 0.2e1 * t129 * t17 + 0.2e1 * t2 * t41 + 0.2e1 * t5 * t29 + 0.2e1 * t3 * t39 + 0.2e1 * t4 * t40 + 0.2e1 * m(5) * (t129 * t3 + t19 * t4) + 0.2e1 * m(6) * (t1 * t13 + t14 * t2 + t21 * t5) + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t124 + (mrSges(4,1) * t124 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * Ifges(4,4) * t65 + 0.2e1 * t68 * t30 + (-t108 - t118) * t64 + t74 * t66) * qJD(3) + t83) * t65 + (mrSges(4,2) * t124 - 0.2e1 * t68 * t11 + (t8 + t9) * t66 + (t6 - t7) * t64 + ((t108 - t116) * t66 + t74 * t64) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) + (-Ifges(5,6) * t64 - 0.2e1 * Ifges(4,4) - t134) * t67 + (t68 ^ 2 * t131 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + t133) * t65) * qJD(3)) * t67; (-m(6) * t5 - t10 - t11 + (t107 * t66 + t106 * t64 + m(6) * (t13 * t66 + t14 * t64) + m(5) * (t129 * t66 - t19 * t64)) * qJD(3)) * t67 + (t90 * t131 + (m(6) * t21 + t29 + t30) * qJD(3) + t126 * t66 + t125 * t64) * t65; 0.4e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t103) * t65 * t101; t5 * t46 + t26 * t29 + t21 * t33 + t44 * t10 - pkin(3) * t11 + m(6) * (t21 * t26 + t44 * t5) + (t61 / 0.2e1 + t58 / 0.2e1 + t60 / 0.2e1 + (-Ifges(4,5) + (-m(5) * pkin(3) + t105) * t68) * qJD(3)) * t65 + (t8 / 0.2e1 + t9 / 0.2e1 + t2 * mrSges(6,2) - t4 * mrSges(5,3) - t86 * t102 + (-t116 / 0.2e1 + t22 / 0.2e1 - t23 / 0.2e1 - t13 * mrSges(6,2) - t129 * mrSges(5,3)) * qJD(4) + t125 * pkin(7)) * t64 + (t3 * mrSges(5,3) + t1 * mrSges(6,2) - t6 / 0.2e1 + t7 / 0.2e1 + t85 * t102 + (t24 / 0.2e1 + t25 / 0.2e1 + t14 * mrSges(6,2) - t19 * mrSges(5,3)) * qJD(4) + t126 * pkin(7)) * t66 + (-t68 * t34 + (t37 / 0.2e1 + t38 / 0.2e1) * t66 + (t35 / 0.2e1 - t36 / 0.2e1) * t64 + (-t68 * mrSges(4,2) - Ifges(4,6) + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t66 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t64) * qJD(3) + (t85 * t64 + t86 * t66) * qJD(4)) * t67; (t46 + t105) * t102 + m(6) * (t44 * t102 + t104) + m(5) * (-pkin(3) * t102 + t104) + (-m(6) * t26 - t34 - t33 + (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t103) * qJD(3)) * t67; -0.2e1 * pkin(3) * t34 + 0.2e1 * t33 * t44 + (-t35 + t36) * t66 + (t37 + t38) * t64 + 0.2e1 * (m(6) * t44 + t46) * t26 + ((t50 + t51) * t66 + (t48 - t49) * t64) * qJD(4); m(6) * (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t13) + t1 * mrSges(6,3) + qJD(5) * t42 + qJ(5) * t18 - t3 * mrSges(5,2) + t4 * mrSges(5,1) - t2 * mrSges(6,1) - pkin(4) * t16 + (-Ifges(5,6) * t66 - t130 * t64) * t98 + t134 * t102 + t83; (-m(6) * t75 - t81 - t82) * t101 + t69 * t65; -t128 * mrSges(6,2) - Ifges(5,6) * t100 + t69 * pkin(7) + t58 + t60 + t61; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t2 + t16; (t64 * t101 + t65 * t99) * m(6); (m(6) * pkin(7) + mrSges(6,2)) * t99; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:01
% EndTime: 2019-12-31 21:11:04
% DurationCPUTime: 1.01s
% Computational Cost: add. (1051->182), mult. (2364->253), div. (0->0), fcn. (1657->6), ass. (0->91)
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t141 = t79 ^ 2 + t82 ^ 2;
t140 = 2 * Ifges(4,4) - 2 * Ifges(5,5);
t110 = pkin(1) * qJD(2);
t83 = cos(qJ(2));
t100 = t83 * t110;
t139 = t141 * t100;
t138 = Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3);
t137 = t82 * pkin(3) + t79 * qJ(4);
t108 = qJD(3) * t82;
t128 = pkin(7) - pkin(8);
t58 = t128 * t79;
t59 = t128 * t82;
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t21 = t58 * t81 - t59 * t78;
t109 = qJD(3) * t79;
t70 = pkin(8) * t109;
t50 = -pkin(7) * t109 + t70;
t51 = qJD(3) * t59;
t8 = qJD(5) * t21 + t50 * t81 + t51 * t78;
t22 = t58 * t78 + t59 * t81;
t9 = -qJD(5) * t22 - t50 * t78 + t51 * t81;
t136 = t9 * mrSges(6,1) - t8 * mrSges(6,2);
t80 = sin(qJ(2));
t63 = pkin(1) * t80 + pkin(7);
t123 = -pkin(8) + t63;
t41 = t123 * t79;
t42 = t123 * t82;
t11 = t41 * t81 - t42 * t78;
t23 = t100 * t82 - t109 * t63 + t70;
t99 = t79 * t100;
t24 = qJD(3) * t42 + t99;
t1 = qJD(5) * t11 + t23 * t81 + t24 * t78;
t12 = t41 * t78 + t42 * t81;
t2 = -qJD(5) * t12 - t23 * t78 + t24 * t81;
t135 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t134 = qJD(3) - qJD(5);
t133 = (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t141) * t100;
t132 = 2 * m(6);
t92 = t78 * t82 - t79 * t81;
t17 = t134 * t92;
t91 = t78 * t79 + t81 * t82;
t18 = t134 * t91;
t5 = mrSges(6,1) * t17 + mrSges(6,2) * t18;
t131 = 0.2e1 * t5;
t19 = mrSges(6,1) * t91 - mrSges(6,2) * t92;
t130 = 0.2e1 * t19;
t94 = t79 * mrSges(5,1) - t82 * mrSges(5,3);
t47 = t94 * qJD(3);
t129 = 0.2e1 * t47;
t118 = t17 * mrSges(6,3);
t84 = -pkin(3) - pkin(4);
t52 = -t78 * qJ(4) + t81 * t84;
t28 = t81 * qJD(4) + qJD(5) * t52;
t117 = t28 * mrSges(6,2);
t53 = t81 * qJ(4) + t78 * t84;
t29 = -t78 * qJD(4) - qJD(5) * t53;
t116 = t29 * mrSges(6,1);
t115 = t91 * mrSges(6,3);
t114 = Ifges(6,5) * t18 - Ifges(6,6) * t17;
t113 = t139 * t63;
t112 = t139 * pkin(7);
t107 = qJD(4) * t82;
t106 = qJD(5) * t78;
t105 = qJD(5) * t81;
t104 = 0.2e1 * mrSges(6,3);
t103 = pkin(2) + t137;
t32 = pkin(3) * t109 - qJ(4) * t108 - t79 * qJD(4);
t102 = m(5) * t107;
t101 = t80 * t110;
t64 = -pkin(1) * t83 - pkin(2);
t56 = -t82 * mrSges(4,1) + t79 * mrSges(4,2);
t95 = t79 * mrSges(4,1) + t82 * mrSges(4,2);
t55 = -t82 * mrSges(5,1) - t79 * mrSges(5,3);
t36 = t64 - t137;
t90 = (-mrSges(3,1) + t56) * t101;
t25 = -pkin(4) * t109 - t32;
t89 = mrSges(5,2) * t108 - t105 * t115 - t78 * t118 + (-t106 * t92 - t18 * t81) * mrSges(6,3);
t88 = -mrSges(5,2) * t137 - Ifges(4,6) * t79;
t87 = -t53 * t118 - t28 * t115 + mrSges(5,2) * t107 + Ifges(5,6) * t109 + (-t18 * t52 + t29 * t92) * mrSges(6,3) - t114 + (Ifges(4,5) + Ifges(5,4)) * t108;
t86 = -0.2e1 * t92 * t18 * Ifges(6,1) + 0.2e1 * t91 * Ifges(6,2) * t17 + (t138 * t82 - t140 * t79) * t109 + (t138 * t79 + t140 * t82) * t108 + 0.2e1 * (t17 * t92 - t18 * t91) * Ifges(6,4);
t85 = -m(5) * t137 + t55 + t56;
t74 = t82 * pkin(4);
t48 = t95 * qJD(3);
t37 = t74 + t103;
t27 = -t36 + t74;
t26 = t101 + t32;
t20 = t25 - t101;
t3 = [0.2e1 * m(4) * (t101 * t64 + t113) + (t1 * t12 + t11 * t2 + t20 * t27) * t132 + 0.2e1 * m(5) * (t26 * t36 + t113) + t86 + 0.2e1 * t64 * t48 + 0.2e1 * t26 * t55 + t36 * t129 + t20 * t130 + t27 * t131 + 0.2e1 * t90 + 0.2e1 * t133 + (-t1 * t91 - t11 * t18 - t12 * t17 + t2 * t92) * t104; (t32 + t26) * t55 + (t37 + t27) * t5 + (t64 - pkin(2)) * t48 + (-t103 + t36) * t47 + (t20 + t25) * t19 + t86 + (-(-t2 - t9) * t92 - (t1 + t8) * t91 + (-t11 - t21) * t18 - (t12 + t22) * t17) * mrSges(6,3) + t133 + t90 + m(4) * (-pkin(2) * t101 + t112) + m(5) * (-t103 * t26 + t32 * t36 + t112) + m(6) * (t1 * t22 + t11 * t9 + t12 * t8 + t2 * t21 + t20 * t37 + t25 * t27); t86 - 0.2e1 * pkin(2) * t48 - t103 * t129 + t37 * t131 + t25 * t130 + (t21 * t9 + t22 * t8 + t25 * t37) * t132 + 0.2e1 * (-m(5) * t103 + t55) * t32 + (-t17 * t22 - t18 * t21 - t8 * t91 + t9 * t92) * t104; t63 * t102 + t87 + m(6) * (t1 * t53 + t11 * t29 + t12 * t28 + t2 * t52) + (t63 * t85 + t88) * qJD(3) + (m(5) * (-pkin(3) * t79 + qJ(4) * t82) - t94 - t95) * t100 - t135; m(6) * (t21 * t29 + t22 * t28 + t52 * t9 + t53 * t8) + t88 * qJD(3) + (qJD(3) * t85 + t102) * pkin(7) + t87 - t136; (t28 * t53 + t29 * t52) * t132 + 0.2e1 * t117 - 0.2e1 * t116 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); m(6) * (t1 * t78 + t2 * t81 + (-t11 * t78 + t12 * t81) * qJD(5)) + (t108 * t63 + t99) * m(5) + t89; m(6) * (t78 * t8 + t81 * t9 + (-t21 * t78 + t22 * t81) * qJD(5)) + m(5) * pkin(7) * t108 + t89; m(6) * (t28 * t78 + t29 * t81 + (-t52 * t78 + t53 * t81) * qJD(5)) + mrSges(6,2) * t105 + mrSges(6,1) * t106; 0; t114 + t135; t114 + t136; t116 - t117; (-t78 * mrSges(6,1) - t81 * mrSges(6,2)) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;

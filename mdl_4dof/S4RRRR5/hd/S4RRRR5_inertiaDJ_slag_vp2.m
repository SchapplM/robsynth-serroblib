% Calculate time derivative of joint inertia matrix for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:23
% EndTime: 2019-12-31 17:27:26
% DurationCPUTime: 1.26s
% Computational Cost: add. (1284->233), mult. (3224->374), div. (0->0), fcn. (2514->6), ass. (0->115)
t84 = sin(qJ(2));
t82 = sin(qJ(4));
t83 = sin(qJ(3));
t85 = cos(qJ(4));
t86 = cos(qJ(3));
t92 = t82 * t83 - t85 * t86;
t52 = t92 * t84;
t87 = cos(qJ(2));
t110 = qJD(2) * t87;
t101 = t86 * t110;
t111 = qJD(2) * t84;
t131 = -Ifges(4,5) * t101 - Ifges(4,3) * t111;
t70 = -pkin(2) * t87 - t84 * pkin(6) - pkin(1);
t115 = t86 * t87;
t77 = pkin(5) * t115;
t48 = t83 * t70 + t77;
t130 = qJD(3) + qJD(4);
t129 = 2 * m(4);
t128 = 2 * m(5);
t127 = -0.2e1 * pkin(1);
t126 = 0.2e1 * pkin(5);
t125 = -t83 / 0.2e1;
t124 = t86 / 0.2e1;
t122 = -pkin(7) - pkin(6);
t121 = pkin(5) * t83;
t120 = Ifges(4,4) * t83;
t119 = Ifges(4,4) * t86;
t118 = Ifges(4,6) * t83;
t117 = t83 * t84;
t116 = t84 * t86;
t34 = t130 * t92;
t60 = t82 * t86 + t83 * t85;
t35 = t130 * t60;
t113 = -Ifges(5,5) * t34 - Ifges(5,6) * t35;
t68 = (pkin(2) * t84 - pkin(6) * t87) * qJD(2);
t112 = t111 * t121 + t86 * t68;
t109 = qJD(3) * t83;
t108 = qJD(3) * t84;
t107 = qJD(3) * t86;
t106 = qJD(4) * t82;
t105 = qJD(4) * t85;
t18 = -t92 * t110 - t35 * t84;
t19 = -t60 * t110 + t130 * t52;
t104 = -Ifges(5,5) * t18 - Ifges(5,6) * t19 - Ifges(5,3) * t111;
t103 = pkin(3) * t109;
t102 = t83 * t108;
t100 = qJD(3) * t122;
t99 = (2 * Ifges(3,4)) + t118;
t98 = -mrSges(4,1) * t86 + mrSges(4,2) * t83;
t97 = mrSges(4,1) * t83 + mrSges(4,2) * t86;
t96 = Ifges(4,1) * t86 - t120;
t72 = Ifges(4,1) * t83 + t119;
t95 = -Ifges(4,2) * t83 + t119;
t71 = Ifges(4,2) * t86 + t120;
t94 = Ifges(4,5) * t83 + Ifges(4,6) * t86;
t23 = t83 * t68 + t70 * t107 + (-t87 * t109 - t86 * t111) * pkin(5);
t24 = -t48 * qJD(3) + t112;
t93 = t23 * t86 - t24 * t83;
t58 = t86 * t70;
t31 = -pkin(7) * t116 + t58 + (-pkin(3) - t121) * t87;
t39 = -pkin(7) * t117 + t48;
t10 = t31 * t85 - t39 * t82;
t11 = t31 * t82 + t39 * t85;
t73 = t122 * t83;
t74 = t122 * t86;
t40 = t73 * t85 + t74 * t82;
t41 = t73 * t82 - t74 * t85;
t66 = t83 * t100;
t67 = t86 * t100;
t21 = t40 * qJD(4) + t66 * t85 + t67 * t82;
t22 = -t41 * qJD(4) - t66 * t82 + t67 * t85;
t91 = t22 * mrSges(5,1) - t21 * mrSges(5,2) + t113;
t88 = t84 * t107 + t83 * t110;
t15 = -t88 * pkin(7) + t23;
t7 = (pkin(3) * t84 - pkin(7) * t115) * qJD(2) + (-t77 + (pkin(7) * t84 - t70) * t83) * qJD(3) + t112;
t2 = t10 * qJD(4) + t15 * t85 + t7 * t82;
t3 = -t11 * qJD(4) - t15 * t82 + t7 * t85;
t90 = t3 * mrSges(5,1) - t2 * mrSges(5,2) - t104;
t89 = t101 - t102;
t81 = Ifges(4,5) * t107;
t78 = -pkin(3) * t86 - pkin(2);
t69 = (pkin(3) * t83 + pkin(5)) * t84;
t65 = -mrSges(4,1) * t87 - mrSges(4,3) * t116;
t64 = mrSges(4,2) * t87 - mrSges(4,3) * t117;
t63 = t96 * qJD(3);
t62 = t95 * qJD(3);
t61 = t97 * qJD(3);
t56 = (-mrSges(5,1) * t82 - mrSges(5,2) * t85) * qJD(4) * pkin(3);
t51 = t60 * t84;
t50 = -Ifges(4,5) * t87 + t96 * t84;
t49 = -Ifges(4,6) * t87 + t95 * t84;
t47 = -t87 * t121 + t58;
t46 = t88 * pkin(3) + pkin(5) * t110;
t45 = -mrSges(4,2) * t111 - t88 * mrSges(4,3);
t44 = mrSges(4,1) * t111 - t89 * mrSges(4,3);
t43 = -mrSges(5,1) * t87 + t52 * mrSges(5,3);
t42 = mrSges(5,2) * t87 - t51 * mrSges(5,3);
t38 = Ifges(5,1) * t60 - Ifges(5,4) * t92;
t37 = Ifges(5,4) * t60 - Ifges(5,2) * t92;
t36 = mrSges(5,1) * t92 + mrSges(5,2) * t60;
t30 = t88 * mrSges(4,1) + t89 * mrSges(4,2);
t29 = mrSges(5,1) * t51 - mrSges(5,2) * t52;
t28 = -t72 * t108 + (Ifges(4,5) * t84 + t96 * t87) * qJD(2);
t27 = -t71 * t108 + (Ifges(4,6) * t84 + t95 * t87) * qJD(2);
t26 = -Ifges(5,1) * t52 - Ifges(5,4) * t51 - Ifges(5,5) * t87;
t25 = -Ifges(5,4) * t52 - Ifges(5,2) * t51 - Ifges(5,6) * t87;
t14 = -Ifges(5,1) * t34 - Ifges(5,4) * t35;
t13 = -Ifges(5,4) * t34 - Ifges(5,2) * t35;
t12 = mrSges(5,1) * t35 - mrSges(5,2) * t34;
t9 = -mrSges(5,2) * t111 + t19 * mrSges(5,3);
t8 = mrSges(5,1) * t111 - t18 * mrSges(5,3);
t6 = -mrSges(5,1) * t19 + mrSges(5,2) * t18;
t5 = Ifges(5,1) * t18 + Ifges(5,4) * t19 + Ifges(5,5) * t111;
t4 = Ifges(5,4) * t18 + Ifges(5,2) * t19 + Ifges(5,6) * t111;
t1 = [0.2e1 * t10 * t8 + 0.2e1 * t11 * t9 + t18 * t26 + t19 * t25 + 0.2e1 * t2 * t42 + 0.2e1 * t23 * t64 + 0.2e1 * t24 * t65 + 0.2e1 * t46 * t29 + 0.2e1 * t3 * t43 - t51 * t4 + 0.2e1 * t47 * t44 + 0.2e1 * t48 * t45 - t52 * t5 + 0.2e1 * t69 * t6 + (t48 * t23 + t47 * t24) * t129 + (t10 * t3 + t11 * t2 + t46 * t69) * t128 + ((mrSges(3,2) * t127 - t83 * t49 + t86 * t50 + t99 * t87) * qJD(2) + t104 + t131) * t87 + (t30 * t126 - t83 * t27 + t86 * t28 + (-t86 * t49 - t83 * t50 + t87 * t94) * qJD(3) + (mrSges(3,1) * t127 - Ifges(5,5) * t52 - Ifges(5,6) * t51 + (Ifges(4,5) * t86 - t99) * t84 + (pkin(5) ^ 2 * t129 + t97 * t126 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3)) * t87) * qJD(2)) * t84; t29 * t103 + (t27 / 0.2e1 + qJD(3) * t50 / 0.2e1) * t86 + t83 * t28 / 0.2e1 + t78 * t6 - (-Ifges(4,6) * t109 + t113 + t81) * t87 / 0.2e1 + (-Ifges(3,6) * qJD(2) + t63 * t124 + t62 * t125 + (-t86 * t71 / 0.2e1 + t72 * t125) * qJD(3) + (qJD(2) * mrSges(3,2) + t61) * pkin(5) + (Ifges(5,5) * t60 - Ifges(5,6) * t92 + t94) * qJD(2) / 0.2e1) * t84 - t92 * t4 / 0.2e1 + (t10 * t34 - t11 * t35 - t2 * t92 - t3 * t60) * mrSges(5,3) + t69 * t12 - t51 * t13 / 0.2e1 - t52 * t14 / 0.2e1 + t60 * t5 / 0.2e1 - t35 * t25 / 0.2e1 + t19 * t37 / 0.2e1 + t18 * t38 / 0.2e1 + t40 * t8 + t41 * t9 + t21 * t42 + t22 * t43 + t46 * t36 - pkin(2) * t30 - t34 * t26 / 0.2e1 + (Ifges(3,5) + t72 * t124 + t71 * t125 + (-m(4) * pkin(2) - mrSges(3,1) + t98) * pkin(5)) * t110 - t49 * t109 / 0.2e1 + (m(4) * (-t47 * t107 - t48 * t109 + t93) + t86 * t45 - t83 * t44 - t65 * t107 - t64 * t109) * pkin(6) + m(5) * (t22 * t10 + t69 * t103 + t21 * t11 + t41 * t2 + t40 * t3 + t46 * t78) + ((-t47 * t86 - t48 * t83) * qJD(3) + t93) * mrSges(4,3); 0.2e1 * t36 * t103 + 0.2e1 * t78 * t12 + (t78 * t103 + t41 * t21 + t40 * t22) * t128 - t35 * t37 - t92 * t13 - t34 * t38 + t60 * t14 - 0.2e1 * pkin(2) * t61 + t83 * t63 - t71 * t109 + (qJD(3) * t72 + t62) * t86 + 0.2e1 * (-t21 * t92 - t22 * t60 + t34 * t40 - t35 * t41) * mrSges(5,3); -Ifges(4,5) * t102 + t24 * mrSges(4,1) - t23 * mrSges(4,2) - t88 * Ifges(4,6) + (m(5) * (-t10 * t106 + t11 * t105 + t2 * t82 + t3 * t85) + t42 * t105 + t82 * t9 - t43 * t106 + t85 * t8) * pkin(3) + t90 - t131; t81 + (t98 * pkin(6) - t118) * qJD(3) + (m(5) * (t21 * t82 + t22 * t85 + (-t40 * t82 + t41 * t85) * qJD(4)) + (t85 * t34 - t82 * t35 + (t60 * t82 - t85 * t92) * qJD(4)) * mrSges(5,3)) * pkin(3) + t91; 0.2e1 * t56; t90; t91; t56; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

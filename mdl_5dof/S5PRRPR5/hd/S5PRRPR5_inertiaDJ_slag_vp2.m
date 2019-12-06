% Calculate time derivative of joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:22
% EndTime: 2019-12-05 16:25:27
% DurationCPUTime: 1.15s
% Computational Cost: add. (1263->215), mult. (3250->344), div. (0->0), fcn. (2981->10), ass. (0->114)
t120 = m(5) * pkin(3);
t61 = sin(pkin(10));
t124 = t61 * t120 - mrSges(5,2);
t101 = cos(pkin(10));
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t51 = -mrSges(6,1) * t67 + mrSges(6,2) * t64;
t86 = t101 * pkin(3);
t57 = -t86 - pkin(4);
t123 = m(6) * t57 - t101 * t120 - mrSges(5,1) + t51;
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t103 = -qJ(4) - pkin(7);
t82 = qJD(3) * t103;
t39 = qJD(4) * t68 + t65 * t82;
t71 = -qJD(4) * t65 + t68 * t82;
t18 = -t101 * t71 + t39 * t61;
t122 = 0.2e1 * t18;
t52 = t103 * t68;
t83 = t101 * t65;
t32 = -t103 * t83 - t52 * t61;
t121 = 0.2e1 * t32;
t119 = t67 / 0.2e1;
t118 = pkin(3) * t61;
t117 = Ifges(6,4) * t64;
t116 = Ifges(6,4) * t67;
t115 = Ifges(6,6) * t64;
t114 = t18 * t32;
t62 = sin(pkin(5));
t66 = sin(qJ(2));
t110 = t62 * t66;
t63 = cos(pkin(5));
t43 = t68 * t110 + t63 * t65;
t69 = cos(qJ(2));
t99 = qJD(2) * t69;
t91 = t62 * t99;
t30 = -t43 * qJD(3) - t65 * t91;
t42 = -t65 * t110 + t63 * t68;
t31 = t42 * qJD(3) + t68 * t91;
t10 = -t101 * t30 + t31 * t61;
t21 = -t101 * t42 + t43 * t61;
t113 = t21 * t10;
t45 = t61 * t68 + t83;
t40 = t45 * qJD(3);
t112 = t40 * mrSges(5,3);
t111 = t61 * t65;
t109 = t62 * t69;
t108 = t64 * mrSges(6,3);
t53 = Ifges(6,2) * t67 + t117;
t107 = t64 * t53;
t106 = t67 * mrSges(6,3);
t72 = t101 * t68 - t111;
t41 = t72 * qJD(3);
t105 = t67 * t41;
t54 = Ifges(6,1) * t64 + t116;
t104 = t67 * t54;
t102 = Ifges(6,5) * t105 + Ifges(6,3) * t40;
t100 = qJD(2) * t66;
t98 = qJD(3) * t65;
t97 = qJD(3) * t68;
t96 = qJD(5) * t64;
t95 = qJD(5) * t67;
t94 = 0.2e1 * t65;
t93 = pkin(3) * t98;
t92 = t62 * t100;
t90 = t45 * t96;
t89 = t45 * t95;
t88 = mrSges(6,3) * t96;
t87 = mrSges(6,3) * t95;
t58 = -pkin(3) * t68 - pkin(2);
t85 = -t96 / 0.2e1;
t23 = t40 * mrSges(5,1) + t41 * mrSges(5,2);
t84 = -(2 * Ifges(5,4)) - t115;
t81 = t62 ^ 2 * t66 * t99;
t80 = mrSges(6,1) * t64 + mrSges(6,2) * t67;
t79 = Ifges(6,1) * t67 - t117;
t78 = -Ifges(6,2) * t64 + t116;
t77 = Ifges(6,5) * t64 + Ifges(6,6) * t67;
t76 = t32 * t10 + t18 * t21;
t25 = -pkin(4) * t72 - pkin(8) * t45 + t58;
t33 = -t101 * t52 + t103 * t111;
t8 = t25 * t67 - t33 * t64;
t9 = t25 * t64 + t33 * t67;
t22 = t101 * t43 + t61 * t42;
t14 = -t67 * t109 - t64 * t22;
t75 = t64 * t109 - t67 * t22;
t74 = t41 * t64 + t89;
t73 = t90 - t105;
t70 = -t30 * t65 + t31 * t68 + (-t42 * t68 - t43 * t65) * qJD(3);
t59 = Ifges(6,5) * t95;
t56 = pkin(8) + t118;
t50 = t79 * qJD(5);
t49 = t78 * qJD(5);
t48 = (mrSges(4,1) * t65 + mrSges(4,2) * t68) * qJD(3);
t47 = t80 * qJD(5);
t29 = -mrSges(5,1) * t72 + mrSges(5,2) * t45;
t27 = -mrSges(6,1) * t72 - t45 * t106;
t26 = mrSges(6,2) * t72 - t45 * t108;
t24 = t80 * t45;
t20 = pkin(4) * t40 - pkin(8) * t41 + t93;
t19 = t101 * t39 + t61 * t71;
t17 = -Ifges(6,5) * t72 + t79 * t45;
t16 = -Ifges(6,6) * t72 + t78 * t45;
t13 = -mrSges(6,2) * t40 - t74 * mrSges(6,3);
t12 = mrSges(6,1) * t40 + t73 * mrSges(6,3);
t11 = t101 * t31 + t61 * t30;
t7 = t74 * mrSges(6,1) - t73 * mrSges(6,2);
t6 = -t73 * Ifges(6,1) - t74 * Ifges(6,4) + t40 * Ifges(6,5);
t5 = -t73 * Ifges(6,4) - t74 * Ifges(6,2) + t40 * Ifges(6,6);
t4 = t75 * qJD(5) - t64 * t11 + t67 * t92;
t3 = t14 * qJD(5) + t67 * t11 + t64 * t92;
t2 = -t9 * qJD(5) - t19 * t64 + t20 * t67;
t1 = t8 * qJD(5) + t19 * t67 + t20 * t64;
t15 = [0.2e1 * m(6) * (t14 * t4 - t3 * t75 + t113) + 0.2e1 * m(5) * (t22 * t11 + t113 - t81) + 0.2e1 * m(4) * (t42 * t30 + t43 * t31 - t81); t10 * t24 + t14 * t12 - t75 * t13 + t21 * t7 + t3 * t26 + t4 * t27 + ((-t23 - t48) * t69 + (-t69 * mrSges(3,2) + (-mrSges(4,1) * t68 + mrSges(4,2) * t65 - mrSges(3,1) + t29) * t66) * qJD(2)) * t62 + m(6) * (-t1 * t75 + t14 * t2 + t3 * t9 + t4 * t8 + t76) + m(5) * (t33 * t11 + t19 * t22 + (t58 * t100 - t69 * t93) * t62 + t76) + (t10 * t45 + t11 * t72 + t21 * t41 - t22 * t40) * mrSges(5,3) + t70 * mrSges(4,3) + (-pkin(2) * t92 + t70 * pkin(7)) * m(4); -0.2e1 * t33 * t112 - 0.2e1 * pkin(2) * t48 + 0.2e1 * t1 * t26 + 0.2e1 * t8 * t12 + 0.2e1 * t9 * t13 + t24 * t122 + 0.2e1 * t2 * t27 + 0.2e1 * t58 * t23 + t7 * t121 + (-Ifges(4,4) * t65 + pkin(3) * t29) * qJD(3) * t94 + 0.2e1 * m(5) * (t19 * t33 + t58 * t93 + t114) + 0.2e1 * m(6) * (t1 * t9 + t2 * t8 + t114) + (mrSges(5,3) * t121 - t64 * t16 + t67 * t17) * t41 + (0.2e1 * Ifges(4,4) * t68 + (Ifges(4,1) - Ifges(4,2)) * t94) * t97 - (-0.2e1 * t19 * mrSges(5,3) + t84 * t41 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t40 + t102) * t72 + (mrSges(5,3) * t122 + 0.2e1 * Ifges(5,1) * t41 - t64 * t5 + t67 * t6 + (Ifges(6,5) * t67 + t84) * t40 + (-t67 * t16 - t64 * t17 + t72 * t77) * qJD(5)) * t45; m(6) * (t3 * t67 - t4 * t64 + (-t14 * t67 + t64 * t75) * qJD(5)) * t56 + t3 * t106 + t75 * t88 - t4 * t108 - t14 * t87 + t21 * t47 - t31 * mrSges(4,2) + t30 * mrSges(4,1) + t124 * t11 + t123 * t10; t64 * t6 / 0.2e1 + t57 * t7 + t32 * t47 - t112 * t118 - t2 * t108 + t17 * t95 / 0.2e1 - t72 * (-Ifges(6,6) * t96 + t59) / 0.2e1 - Ifges(4,6) * t98 - t8 * t87 - t9 * t88 - t53 * t89 / 0.2e1 + Ifges(4,5) * t97 + t1 * t106 + t16 * t85 + t5 * t119 + (-Ifges(5,6) + t77 / 0.2e1) * t40 + t124 * t19 + (-mrSges(4,1) * t97 + mrSges(4,2) * t98) * pkin(7) + (-t64 * t49 / 0.2e1 + t54 * t85 + t50 * t119) * t45 + t123 * t18 + (m(6) * (t1 * t67 - t2 * t64 + (-t64 * t9 - t67 * t8) * qJD(5)) + t67 * t13 - t64 * t12 - t27 * t95 - t26 * t96) * t56 + (Ifges(5,5) + t104 / 0.2e1 - t107 / 0.2e1 - mrSges(5,3) * t86) * t41; 0.2e1 * t47 * t57 + t49 * t67 + t50 * t64 + (t104 - t107) * qJD(5); m(6) * (t3 * t64 + t4 * t67 + (-t14 * t64 - t67 * t75) * qJD(5)) + m(5) * t92; m(6) * (t1 * t64 + t2 * t67 + (-t64 * t8 + t67 * t9) * qJD(5)) + t26 * t95 + t64 * t13 - t27 * t96 + t67 * t12 + m(5) * t93 + t23; 0; 0; mrSges(6,1) * t4 - mrSges(6,2) * t3; mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,5) * t90 - t74 * Ifges(6,6) + t102; t59 + (t51 * t56 - t115) * qJD(5); -t47; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;

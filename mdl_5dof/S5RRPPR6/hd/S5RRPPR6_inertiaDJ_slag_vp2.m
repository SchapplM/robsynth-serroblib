% Calculate time derivative of joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:37
% EndTime: 2019-12-31 19:31:41
% DurationCPUTime: 1.26s
% Computational Cost: add. (2071->231), mult. (4738->371), div. (0->0), fcn. (4395->8), ass. (0->100)
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t91 = sin(qJ(5));
t93 = cos(qJ(5));
t96 = t87 * t91 - t89 * t93;
t69 = t96 * qJD(5);
t121 = 2 * m(5);
t120 = 2 * m(6);
t86 = t89 ^ 2;
t94 = cos(qJ(2));
t84 = -pkin(2) * t94 - pkin(1);
t118 = 0.2e1 * t84;
t88 = sin(pkin(8));
t81 = pkin(2) * t88 + qJ(4);
t117 = pkin(7) + t81;
t116 = Ifges(5,4) * t87;
t115 = Ifges(5,4) * t89;
t92 = sin(qJ(2));
t105 = -qJ(3) - pkin(6);
t98 = qJD(2) * t105;
t66 = qJD(3) * t94 + t92 * t98;
t90 = cos(pkin(8));
t95 = -t92 * qJD(3) + t94 * t98;
t38 = t66 * t88 - t90 * t95;
t79 = t105 * t92;
t80 = t105 * t94;
t57 = -t90 * t79 - t80 * t88;
t114 = t38 * t57;
t75 = t88 * t94 + t90 * t92;
t67 = t75 * qJD(2);
t113 = t67 * Ifges(5,5);
t112 = t67 * Ifges(5,6);
t73 = t88 * t92 - t90 * t94;
t68 = t73 * qJD(2);
t111 = t68 * t87;
t110 = t68 * t89;
t109 = t75 * t87;
t108 = t87 * mrSges(5,3);
t107 = t87 * Ifges(5,2);
t106 = t89 * mrSges(5,3);
t100 = pkin(2) * qJD(2) * t92;
t30 = pkin(3) * t67 + qJ(4) * t68 - qJD(4) * t75 + t100;
t39 = t90 * t66 + t88 * t95;
t13 = t87 * t30 + t89 * t39;
t51 = t73 * pkin(3) - t75 * qJ(4) + t84;
t58 = t79 * t88 - t80 * t90;
t24 = t87 * t51 + t89 * t58;
t37 = -mrSges(5,1) * t111 - mrSges(5,2) * t110;
t76 = t87 * t93 + t89 * t91;
t70 = t76 * qJD(5);
t104 = -Ifges(6,5) * t69 - Ifges(6,6) * t70;
t102 = 0.2e1 * t94;
t21 = t68 * t96 - t70 * t75;
t22 = t68 * t76 + t69 * t75;
t101 = Ifges(6,5) * t21 + Ifges(6,6) * t22 + Ifges(6,3) * t67;
t83 = -pkin(2) * t90 - pkin(3);
t99 = t67 * mrSges(4,1) - t68 * mrSges(4,2);
t7 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t12 = t89 * t30 - t39 * t87;
t23 = t89 * t51 - t58 * t87;
t97 = Ifges(5,5) * t89 - Ifges(5,6) * t87;
t14 = -pkin(7) * t75 * t89 + pkin(4) * t73 + t23;
t20 = -pkin(7) * t109 + t24;
t3 = t14 * t93 - t20 * t91;
t4 = t14 * t91 + t20 * t93;
t71 = t117 * t87;
t72 = t117 * t89;
t46 = -t71 * t93 - t72 * t91;
t47 = -t71 * t91 + t72 * t93;
t78 = -pkin(4) * t89 + t83;
t55 = Ifges(6,1) * t76 - Ifges(6,4) * t96;
t54 = Ifges(6,4) * t76 - Ifges(6,2) * t96;
t53 = mrSges(5,1) * t73 - t106 * t75;
t52 = -mrSges(5,2) * t73 - t108 * t75;
t50 = -Ifges(6,1) * t69 - Ifges(6,4) * t70;
t49 = -Ifges(6,4) * t69 - Ifges(6,2) * t70;
t48 = mrSges(6,1) * t70 - mrSges(6,2) * t69;
t45 = mrSges(5,1) * t67 + t106 * t68;
t44 = -mrSges(5,2) * t67 + t108 * t68;
t41 = t96 * t75;
t40 = t76 * t75;
t36 = pkin(4) * t109 + t57;
t34 = -qJD(4) * t76 - qJD(5) * t47;
t33 = -qJD(4) * t96 + qJD(5) * t46;
t32 = mrSges(6,1) * t73 + mrSges(6,3) * t41;
t31 = -mrSges(6,2) * t73 - mrSges(6,3) * t40;
t29 = t113 - (t89 * Ifges(5,1) - t116) * t68;
t28 = t112 - (-t107 + t115) * t68;
t25 = -pkin(4) * t111 + t38;
t16 = -Ifges(6,1) * t41 - Ifges(6,4) * t40 + Ifges(6,5) * t73;
t15 = -Ifges(6,4) * t41 - Ifges(6,2) * t40 + Ifges(6,6) * t73;
t11 = -mrSges(6,2) * t67 + mrSges(6,3) * t22;
t10 = mrSges(6,1) * t67 - mrSges(6,3) * t21;
t9 = pkin(7) * t111 + t13;
t8 = pkin(4) * t67 + pkin(7) * t110 + t12;
t6 = Ifges(6,1) * t21 + Ifges(6,4) * t22 + t67 * Ifges(6,5);
t5 = Ifges(6,4) * t21 + Ifges(6,2) * t22 + t67 * Ifges(6,6);
t2 = -qJD(5) * t4 + t8 * t93 - t9 * t91;
t1 = qJD(5) * t3 + t8 * t91 + t9 * t93;
t17 = [t99 * t118 + t67 * (-Ifges(6,5) * t41 - Ifges(6,6) * t40) + 0.2e1 * t13 * t52 + 0.2e1 * t12 * t53 + 0.2e1 * t57 * t37 - t40 * t5 + 0.2e1 * t25 * (mrSges(6,1) * t40 - mrSges(6,2) * t41) - t41 * t6 + 0.2e1 * t24 * t44 + 0.2e1 * t23 * t45 + 0.2e1 * t1 * t31 + 0.2e1 * t2 * t32 + 0.2e1 * t36 * t7 + t22 * t15 + t21 * t16 + 0.2e1 * t3 * t10 + 0.2e1 * t4 * t11 + 0.2e1 * (-t57 * t68 - t58 * t67) * mrSges(4,3) + 0.2e1 * m(4) * (t39 * t58 + t114) + (t1 * t4 + t2 * t3 + t25 * t36) * t120 + (t12 * t23 + t13 * t24 + t114) * t121 + (-0.2e1 * t39 * mrSges(4,3) - 0.2e1 * (-Ifges(4,4) + t97) * t68 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(6,3)) * t67 + t101) * t73 + (-t87 * t28 + t89 * t29 + (-0.2e1 * Ifges(4,4) + t97) * t67 - (Ifges(5,1) * t86 + (2 * Ifges(4,1)) + (t107 - 0.2e1 * t115) * t87) * t68 + 0.2e1 * (mrSges(5,1) * t87 + mrSges(5,2) * t89 + mrSges(4,3)) * t38) * t75 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t94) * t102 + (m(4) * pkin(2) * t118 - 0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (mrSges(4,1) * t73 + mrSges(4,2) * t75) - 0.2e1 * Ifges(3,4) * t92 + (Ifges(3,1) - Ifges(3,2)) * t102) * t92) * qJD(2); m(6) * (t1 * t47 + t2 * t46 + t25 * t78 + t3 * t34 + t33 * t4) + t83 * t37 + t76 * t6 / 0.2e1 + t78 * t7 - t69 * t16 / 0.2e1 - t70 * t15 / 0.2e1 - Ifges(4,6) * t67 - t41 * t50 / 0.2e1 + t22 * t54 / 0.2e1 + t21 * t55 / 0.2e1 - t39 * mrSges(4,2) + t46 * t10 + t47 * t11 + t36 * t48 - t40 * t49 / 0.2e1 + t33 * t31 + t34 * t32 - t38 * mrSges(4,1) - (t89 * (Ifges(5,1) * t87 + t115) / 0.2e1 - t87 * (Ifges(5,2) * t89 + t116) / 0.2e1 + Ifges(4,5)) * t68 + (t81 * t44 + qJD(4) * t52 + t13 * mrSges(5,3) + t112 / 0.2e1 - t38 * mrSges(5,1) + t28 / 0.2e1) * t89 + (-qJD(4) * t53 - t12 * mrSges(5,3) - t81 * t45 + t113 / 0.2e1 + t38 * mrSges(5,2) + t29 / 0.2e1) * t87 + t73 * t104 / 0.2e1 + (m(4) * (-t38 * t90 + t39 * t88) + (-t88 * t67 + t90 * t68) * mrSges(4,3)) * pkin(2) - t96 * t5 / 0.2e1 + t67 * (Ifges(6,5) * t76 - Ifges(6,6) * t96) / 0.2e1 + t25 * (mrSges(6,1) * t96 + mrSges(6,2) * t76) + (-t1 * t96 - t2 * t76 + t3 * t69 - t4 * t70) * mrSges(6,3) + m(5) * (t38 * t83 + (-t12 * t87 + t13 * t89) * t81 + (-t23 * t87 + t24 * t89) * qJD(4)) + (Ifges(3,5) * t94 - Ifges(3,6) * t92 + (-mrSges(3,1) * t94 + mrSges(3,2) * t92) * pkin(6)) * qJD(2); 0.2e1 * t78 * t48 + (t33 * t47 + t34 * t46) * t120 - t69 * t55 + t76 * t50 - t70 * t54 - t96 * t49 + 0.2e1 * (-t33 * t96 - t34 * t76 + t46 * t69 - t47 * t70) * mrSges(6,3) + (t121 * t81 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t87 ^ 2 + t86); m(4) * t100 - t96 * t10 + t76 * t11 - t69 * t31 - t70 * t32 + t87 * t44 + t89 * t45 + m(6) * (t1 * t76 - t2 * t96 - t3 * t70 - t4 * t69) + m(5) * (t12 * t89 + t13 * t87) + t99; m(6) * (t33 * t76 - t34 * t96 - t46 * t70 - t47 * t69); (-t69 * t76 + t70 * t96) * t120; m(5) * t38 + m(6) * t25 + t37 + t7; t48; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t101; mrSges(6,1) * t34 - mrSges(6,2) * t33 + t104; -t48; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:35
% EndTime: 2019-12-05 17:09:39
% DurationCPUTime: 1.08s
% Computational Cost: add. (1268->168), mult. (3112->259), div. (0->0), fcn. (2641->8), ass. (0->88)
t124 = qJD(2) + qJD(3);
t72 = sin(qJ(3));
t73 = sin(qJ(2));
t76 = cos(qJ(3));
t77 = cos(qJ(2));
t50 = t72 * t73 - t76 * t77;
t33 = t124 * t50;
t71 = sin(qJ(4));
t75 = cos(qJ(4));
t99 = t71 ^ 2 + t75 ^ 2;
t85 = t99 * t33;
t58 = -mrSges(5,1) * t75 + mrSges(5,2) * t71;
t130 = -mrSges(4,1) + t58;
t129 = Ifges(5,1) - Ifges(5,2);
t52 = t72 * t77 + t76 * t73;
t70 = sin(qJ(5));
t74 = cos(qJ(5));
t81 = t70 * t71 - t74 * t75;
t22 = t81 * t52;
t63 = pkin(2) * t72 + pkin(7);
t114 = -pkin(8) - t63;
t47 = t114 * t71;
t67 = t75 * pkin(8);
t48 = t63 * t75 + t67;
t23 = t47 * t74 - t48 * t70;
t84 = qJD(4) * t114;
t98 = pkin(2) * qJD(3);
t89 = t76 * t98;
t39 = t71 * t84 + t75 * t89;
t40 = -t71 * t89 + t75 * t84;
t6 = qJD(5) * t23 + t39 * t74 + t40 * t70;
t24 = t47 * t70 + t48 * t74;
t7 = -qJD(5) * t24 - t39 * t70 + t40 * t74;
t127 = t7 * mrSges(6,1) - t6 * mrSges(6,2);
t117 = -pkin(8) - pkin(7);
t60 = t117 * t71;
t61 = pkin(7) * t75 + t67;
t36 = t60 * t74 - t61 * t70;
t86 = qJD(4) * t117;
t54 = t71 * t86;
t55 = t75 * t86;
t17 = qJD(5) * t36 + t54 * t74 + t55 * t70;
t37 = t60 * t70 + t61 * t74;
t18 = -qJD(5) * t37 - t54 * t70 + t55 * t74;
t126 = t18 * mrSges(6,1) - t17 * mrSges(6,2);
t125 = t99 * t76;
t123 = qJD(4) + qJD(5);
t122 = (mrSges(5,3) * t99 - mrSges(4,2)) * t76 + t130 * t72;
t121 = 2 * m(6);
t31 = t123 * t81;
t51 = t70 * t75 + t71 * t74;
t32 = t123 * t51;
t12 = mrSges(6,1) * t32 - mrSges(6,2) * t31;
t120 = 0.2e1 * t12;
t35 = mrSges(6,1) * t81 + mrSges(6,2) * t51;
t119 = 0.2e1 * t35;
t118 = m(4) / 0.2e1;
t116 = pkin(2) * t76;
t111 = Ifges(5,6) * t71;
t109 = t32 * mrSges(6,3);
t34 = t124 * t52;
t19 = t50 * t34;
t106 = t50 * t72;
t100 = -Ifges(6,5) * t31 - Ifges(6,6) * t32;
t97 = pkin(4) * qJD(5);
t96 = qJD(4) * t71;
t95 = qJD(4) * t75;
t94 = qJD(5) * t70;
t93 = qJD(5) * t74;
t92 = 0.2e1 * mrSges(6,3);
t91 = t74 * t31 * mrSges(6,3);
t90 = mrSges(6,3) * t97;
t88 = pkin(4) * t96;
t3 = -t32 * t52 + t81 * t33;
t4 = t123 * t22 + t51 * t33;
t87 = t4 * mrSges(6,1) - t3 * mrSges(6,2);
t65 = -pkin(4) * t75 - pkin(3);
t82 = mrSges(5,1) * t71 + mrSges(5,2) * t75;
t80 = -t74 * t81 * t90 + Ifges(5,5) * t95 + t100 + (-pkin(4) * t109 + t51 * t90) * t70;
t79 = -0.2e1 * t51 * t31 * Ifges(6,1) + 0.2e1 * t32 * Ifges(6,2) * t81 + (-0.2e1 * Ifges(5,4) * t71 + t129 * t75) * t96 + (0.2e1 * Ifges(5,4) * t75 + t129 * t71) * t95 + 0.2e1 * (t31 * t81 - t51 * t32) * Ifges(6,4);
t21 = t51 * t52;
t53 = t82 * qJD(4);
t78 = t33 * mrSges(4,2) + t22 * t109 + (t53 + t12) * t50 + (-t21 * t31 - t3 * t81 - t4 * t51) * mrSges(6,3) - t85 * mrSges(5,3) + (t35 + t130) * t34;
t64 = -pkin(3) - t116;
t57 = t65 - t116;
t56 = t72 * t98 + t88;
t46 = (-mrSges(6,1) * t70 - mrSges(6,2) * t74) * t97;
t1 = [0.2e1 * m(6) * (-t21 * t4 - t22 * t3 + t19) + 0.2e1 * m(5) * (-t52 * t85 + t19) + 0.2e1 * m(4) * (-t33 * t52 + t19); 0.2e1 * ((-t33 * t72 - t34 * t76) * t118 + (m(5) * (t125 * t52 + t106) / 0.2e1 + (t52 * t76 + t106) * t118) * qJD(3)) * pkin(2) + m(5) * (t34 * t64 - t63 * t85) + m(6) * (-t21 * t7 - t22 * t6 + t23 * t4 + t24 * t3 + t34 * t57 + t50 * t56) + (-mrSges(3,1) * t73 - mrSges(3,2) * t77) * qJD(2) + t78; (t23 * t7 + t24 * t6 + t56 * t57) * t121 + t56 * t119 + t57 * t120 + 0.2e1 * t64 * t53 + (t23 * t31 - t24 * t32 - t51 * t7 - t6 * t81) * t92 + 0.2e1 * (m(5) * (t125 * t63 + t64 * t72) + t122) * t98 + t79; m(6) * (-t17 * t22 - t18 * t21 + t3 * t37 + t34 * t65 + t36 * t4 + t50 * t88) + t78 + m(5) * (-pkin(3) * t34 - pkin(7) * t85); m(6) * (t17 * t24 + t18 * t23 + t36 * t7 + t37 * t6 + t56 * t65 + t57 * t88) + (-pkin(3) + t64) * t53 + (t56 + t88) * t35 + (t57 + t65) * t12 + (m(5) * (-pkin(3) * t72 + t125 * pkin(7)) + t122) * t98 + ((-t18 - t7) * t51 - (t17 + t6) * t81 - (t24 + t37) * t32 - (-t23 - t36) * t31) * mrSges(6,3) + t79; t88 * t119 + t65 * t120 - 0.2e1 * pkin(3) * t53 + (t17 * t37 + t18 * t36 + t65 * t88) * t121 + (-t17 * t81 - t18 * t51 + t31 * t36 - t32 * t37) * t92 + t79; (t33 * t75 + t52 * t96) * mrSges(5,2) + (t33 * t71 - t52 * t95) * mrSges(5,1) + m(6) * (t3 * t70 + t4 * t74 + (t21 * t70 - t22 * t74) * qJD(5)) * pkin(4) + t87; -t82 * t89 + (t58 * t63 - t111) * qJD(4) + (t91 + m(6) * (-t23 * t94 + t24 * t93 + t6 * t70 + t7 * t74)) * pkin(4) + t80 + t127; (t58 * pkin(7) - t111) * qJD(4) + (t91 + m(6) * (t17 * t70 + t18 * t74 - t36 * t94 + t37 * t93)) * pkin(4) + t80 + t126; 0.2e1 * t46; t87; t100 + t127; t100 + t126; t46; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

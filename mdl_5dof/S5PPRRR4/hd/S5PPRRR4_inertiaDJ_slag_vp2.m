% Calculate time derivative of joint inertia matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:39
% EndTime: 2019-12-05 15:18:44
% DurationCPUTime: 1.22s
% Computational Cost: add. (1074->209), mult. (3378->347), div. (0->0), fcn. (3329->12), ass. (0->109)
t52 = sin(pkin(11));
t54 = sin(pkin(5));
t57 = cos(pkin(5));
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t55 = cos(pkin(11));
t56 = cos(pkin(6));
t96 = t55 * t56;
t53 = sin(pkin(6));
t97 = t53 * t63;
t123 = (-t52 * t60 + t63 * t96) * t54 + t57 * t97;
t121 = m(6) * pkin(9) + mrSges(6,3);
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t43 = -mrSges(6,1) * t61 + mrSges(6,2) * t58;
t120 = -m(6) * pkin(4) - mrSges(5,1) + t43;
t119 = 0.2e1 * m(6);
t118 = 2 * pkin(8);
t117 = m(5) / 0.2e1;
t116 = m(6) / 0.2e1;
t115 = m(5) * pkin(3);
t98 = t53 * t60;
t19 = t57 * t98 + (t52 * t63 + t60 * t96) * t54;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t70 = -t53 * t54 * t55 + t56 * t57;
t10 = t19 * t62 + t70 * t59;
t15 = t123 * qJD(3);
t3 = t10 * qJD(4) + t15 * t59;
t9 = t19 * t59 - t70 * t62;
t114 = t3 * t9;
t107 = Ifges(6,4) * t58;
t45 = Ifges(6,2) * t61 + t107;
t113 = -t45 / 0.2e1;
t112 = -t58 / 0.2e1;
t111 = pkin(8) * t62;
t110 = t3 * t59;
t4 = -t9 * qJD(4) + t15 * t62;
t109 = t4 * t62;
t108 = mrSges(6,3) * t59;
t106 = Ifges(6,4) * t61;
t105 = Ifges(6,6) * t58;
t16 = qJD(3) * t19;
t104 = t16 * t123;
t31 = -t62 * t56 + t59 * t98;
t91 = qJD(3) * t63;
t82 = t53 * t91;
t20 = -t31 * qJD(4) + t62 * t82;
t103 = t20 * t62;
t32 = t56 * t59 + t62 * t98;
t21 = t32 * qJD(4) + t59 * t82;
t102 = t21 * t59;
t101 = t31 * t21;
t42 = -pkin(4) * t62 - pkin(9) * t59 - pkin(3);
t100 = t42 * t58;
t95 = t62 * Ifges(6,6);
t94 = -mrSges(5,1) * t62 + mrSges(5,2) * t59 - mrSges(4,1);
t89 = qJD(4) * t62;
t80 = t61 * t89;
t90 = qJD(4) * t59;
t93 = -Ifges(6,5) * t80 - Ifges(6,3) * t90;
t92 = qJD(3) * t60;
t88 = qJD(5) * t59;
t87 = qJD(5) * t61;
t86 = qJD(5) * t62;
t83 = t53 * t92;
t81 = t58 * t88;
t79 = (2 * Ifges(5,4)) + t105;
t41 = (pkin(4) * t59 - pkin(9) * t62) * qJD(4);
t11 = t42 * t87 + t41 * t58 + (-t58 * t86 - t61 * t90) * pkin(8);
t26 = -t58 * t111 + t42 * t61;
t78 = -t26 * qJD(5) + t11;
t77 = t21 * t9 + t31 * t3;
t76 = mrSges(6,1) * t58 + mrSges(6,2) * t61;
t75 = Ifges(6,1) * t61 - t107;
t46 = Ifges(6,1) * t58 + t106;
t74 = -Ifges(6,2) * t58 + t106;
t73 = Ifges(6,5) * t58 + Ifges(6,6) * t61;
t6 = t10 * t61 - t123 * t58;
t5 = -t10 * t58 - t123 * t61;
t72 = t9 * t89 + t110;
t22 = -t58 * t32 - t61 * t97;
t71 = -t61 * t32 + t58 * t97;
t69 = -t123 * t92 - t16 * t63;
t68 = t31 * t89 + t102;
t67 = t80 - t81;
t66 = t58 * t89 + t59 * t87;
t51 = Ifges(6,5) * t87;
t40 = -mrSges(6,1) * t62 - t61 * t108;
t39 = mrSges(6,2) * t62 - t58 * t108;
t37 = t75 * qJD(5);
t36 = t74 * qJD(5);
t35 = (mrSges(5,1) * t59 + mrSges(5,2) * t62) * qJD(4);
t34 = t76 * qJD(5);
t33 = t76 * t59;
t30 = -Ifges(6,5) * t62 + t75 * t59;
t29 = t74 * t59 - t95;
t27 = t61 * t111 + t100;
t25 = -mrSges(6,2) * t90 - t66 * mrSges(6,3);
t24 = mrSges(6,1) * t90 - t67 * mrSges(6,3);
t17 = t66 * mrSges(6,1) + t67 * mrSges(6,2);
t14 = -t46 * t88 + (Ifges(6,5) * t59 + t75 * t62) * qJD(4);
t13 = -t45 * t88 + (Ifges(6,6) * t59 + t74 * t62) * qJD(4);
t12 = -qJD(5) * t100 + t41 * t61 + (t58 * t90 - t61 * t86) * pkin(8);
t8 = t71 * qJD(5) - t58 * t20 + t61 * t83;
t7 = t22 * qJD(5) + t61 * t20 + t58 * t83;
t2 = t5 * qJD(5) + t16 * t58 + t4 * t61;
t1 = -t6 * qJD(5) + t16 * t61 - t4 * t58;
t18 = [0.2e1 * m(6) * (t1 * t5 + t2 * t6 + t114) + 0.2e1 * m(5) * (t10 * t4 - t104 + t114) + 0.2e1 * m(4) * (t15 * t19 - t104); m(6) * (t1 * t22 - t2 * t71 + t5 * t8 + t6 * t7 + t77) + m(5) * (t20 * t10 + t32 * t4 + t77) + 0.2e1 * (t69 * t117 + m(4) * (t15 * t60 + t19 * t91 + t69) / 0.2e1) * t53; 0.2e1 * m(5) * (-t53 ^ 2 * t60 * t91 + t32 * t20 + t101) + 0.2e1 * m(6) * (t22 * t8 - t7 * t71 + t101); -t15 * mrSges(4,2) + t1 * t40 + t9 * t17 - t123 * t35 + t2 * t39 + t5 * t24 + t6 * t25 + t3 * t33 + m(6) * (t1 * t26 + t11 * t6 + t12 * t5 + t2 * t27) + (t72 * t116 + (-t10 * t90 + t109 + t72) * t117) * t118 + (t110 + t109 + (-t10 * t59 + t62 * t9) * qJD(4)) * mrSges(5,3) + (t94 - t115) * t16; t31 * t17 + t21 * t33 + t22 * t24 - t71 * t25 + t7 * t39 + t8 * t40 + (-t63 * t35 + (-mrSges(4,2) * t63 + t94 * t60) * qJD(3)) * t53 - t83 * t115 + m(6) * (-t11 * t71 + t12 * t22 + t26 * t8 + t27 * t7) + ((-t32 * t90 + t103 + t68) * t117 + t68 * t116) * t118 + (t103 + t102 + (t31 * t62 - t32 * t59) * qJD(4)) * mrSges(5,3); 0.2e1 * t12 * t40 + 0.2e1 * t26 * t24 + (t11 * t27 + t12 * t26) * t119 + 0.2e1 * t11 * t39 + 0.2e1 * t27 * t25 - 0.2e1 * pkin(3) * t35 + ((t33 * t118 - t58 * t29 + t61 * t30 + t79 * t62) * qJD(4) + t93) * t62 + (t17 * t118 - t58 * t13 + t61 * t14 + (-t61 * t29 - t58 * t30 + t62 * t73) * qJD(5) + ((Ifges(6,5) * t61 - t79) * t59 + ((pkin(8) ^ 2) * t119 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) - Ifges(6,3)) * t62) * qJD(4)) * t59; -t4 * mrSges(5,2) + t9 * t34 + t121 * (-t1 * t58 + t2 * t61 + (-t5 * t61 - t58 * t6) * qJD(5)) + t120 * t3; -t20 * mrSges(5,2) + t31 * t34 + t121 * (-t8 * t58 + t7 * t61 + (-t22 * t61 + t58 * t71) * qJD(5)) + t120 * t21; -pkin(4) * t17 + (-t51 / 0.2e1 + (t120 * pkin(8) + Ifges(5,5)) * qJD(4)) * t62 + (t14 / 0.2e1 - t12 * mrSges(6,3) + t89 * t113 + (t95 / 0.2e1 - t29 / 0.2e1 - t27 * mrSges(6,3)) * qJD(5) + (m(6) * (-t27 * qJD(5) - t12) - t24 - qJD(5) * t39) * pkin(9)) * t58 + (qJD(5) * t30 / 0.2e1 + t13 / 0.2e1 + t46 * t89 / 0.2e1 + t78 * mrSges(6,3) + (m(6) * t78 - qJD(5) * t40 + t25) * pkin(9)) * t61 + (t36 * t112 + t61 * t37 / 0.2e1 + (t46 * t112 + t61 * t113) * qJD(5) + pkin(8) * t34 + (-Ifges(5,6) + t73 / 0.2e1 + pkin(8) * mrSges(5,2)) * qJD(4)) * t59; -0.2e1 * pkin(4) * t34 + t36 * t61 + t37 * t58 + (-t45 * t58 + t46 * t61) * qJD(5); mrSges(6,1) * t1 - mrSges(6,2) * t2; mrSges(6,1) * t8 - mrSges(6,2) * t7; mrSges(6,1) * t12 - mrSges(6,2) * t11 - Ifges(6,5) * t81 - t66 * Ifges(6,6) - t93; t51 + (t43 * pkin(9) - t105) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t18(1), t18(2), t18(4), t18(7), t18(11); t18(2), t18(3), t18(5), t18(8), t18(12); t18(4), t18(5), t18(6), t18(9), t18(13); t18(7), t18(8), t18(9), t18(10), t18(14); t18(11), t18(12), t18(13), t18(14), t18(15);];
Mq = res;

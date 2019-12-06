% Calculate time derivative of joint inertia matrix for
% S5RRRPR1
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:49
% EndTime: 2019-12-05 18:37:53
% DurationCPUTime: 1.30s
% Computational Cost: add. (3497->189), mult. (7784->308), div. (0->0), fcn. (7264->8), ass. (0->87)
t109 = -pkin(7) - pkin(6);
t89 = sin(qJ(2));
t79 = t109 * t89;
t92 = cos(qJ(2));
t80 = t109 * t92;
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t54 = t91 * t79 + t80 * t88;
t55 = t88 * t79 - t91 * t80;
t113 = qJD(2) + qJD(3);
t112 = 2 * m(5);
t111 = 2 * m(6);
t83 = -pkin(2) * t92 - pkin(1);
t110 = 0.2e1 * t83;
t85 = sin(pkin(9));
t108 = pkin(3) * t85;
t103 = t85 * t88;
t82 = pkin(2) * t91 + pkin(3);
t86 = cos(pkin(9));
t63 = -pkin(2) * t103 + t86 * t82;
t60 = pkin(4) + t63;
t102 = t86 * t88;
t65 = pkin(2) * t102 + t82 * t85;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t43 = t60 * t90 - t65 * t87;
t101 = pkin(2) * qJD(3);
t61 = (-t85 * t91 - t102) * t101;
t62 = (t86 * t91 - t103) * t101;
t27 = t43 * qJD(5) + t61 * t87 + t62 * t90;
t107 = t27 * mrSges(6,2);
t106 = t61 * mrSges(5,1);
t105 = t62 * mrSges(5,2);
t98 = qJD(2) * t109;
t76 = t89 * t98;
t77 = t92 * t98;
t37 = t54 * qJD(3) + t91 * t76 + t88 * t77;
t72 = t88 * t92 + t91 * t89;
t53 = t113 * t72;
t71 = -t88 * t89 + t91 * t92;
t22 = -t53 * qJ(4) + qJD(4) * t71 + t37;
t38 = -qJD(3) * t55 - t76 * t88 + t91 * t77;
t52 = t113 * t71;
t23 = -t52 * qJ(4) - qJD(4) * t72 + t38;
t15 = t86 * t22 + t85 * t23;
t45 = -qJ(4) * t72 + t54;
t46 = qJ(4) * t71 + t55;
t26 = t85 * t45 + t86 * t46;
t100 = 0.2e1 * t92;
t48 = t71 * t86 - t72 * t85;
t49 = t71 * t85 + t72 * t86;
t29 = t48 * t90 - t49 * t87;
t34 = -t52 * t85 - t53 * t86;
t35 = t52 * t86 - t53 * t85;
t12 = t29 * qJD(5) + t34 * t87 + t35 * t90;
t30 = t48 * t87 + t49 * t90;
t13 = -t30 * qJD(5) + t34 * t90 - t35 * t87;
t99 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t47 = qJD(2) * t89 * pkin(2) + t53 * pkin(3);
t97 = -t34 * mrSges(5,1) + t35 * mrSges(5,2);
t81 = pkin(3) * t86 + pkin(4);
t66 = t90 * t108 + t81 * t87;
t58 = t66 * qJD(5);
t56 = t58 * mrSges(6,1);
t64 = -t87 * t108 + t81 * t90;
t57 = t64 * qJD(5);
t96 = -t57 * mrSges(6,2) - t56;
t14 = -t22 * t85 + t86 * t23;
t25 = t86 * t45 - t46 * t85;
t4 = -t35 * pkin(8) + t14;
t5 = pkin(8) * t34 + t15;
t16 = -t49 * pkin(8) + t25;
t17 = pkin(8) * t48 + t26;
t6 = t16 * t90 - t17 * t87;
t2 = t6 * qJD(5) + t4 * t87 + t5 * t90;
t7 = t16 * t87 + t17 * t90;
t3 = -t7 * qJD(5) + t4 * t90 - t5 * t87;
t95 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t12 + Ifges(6,6) * t13;
t44 = t60 * t87 + t65 * t90;
t59 = -t71 * pkin(3) + t83;
t94 = (-mrSges(4,1) * t88 - mrSges(4,2) * t91) * t101;
t93 = t38 * mrSges(4,1) + t14 * mrSges(5,1) - t37 * mrSges(4,2) - t15 * mrSges(5,2) + Ifges(4,5) * t52 + Ifges(5,5) * t35 - Ifges(4,6) * t53 + Ifges(5,6) * t34 + t95;
t39 = -t48 * pkin(4) + t59;
t28 = -qJD(5) * t44 + t61 * t90 - t62 * t87;
t24 = t28 * mrSges(6,1);
t21 = -t34 * pkin(4) + t47;
t1 = [0.2e1 * t29 * Ifges(6,2) * t13 + 0.2e1 * t30 * t12 * Ifges(6,1) + 0.2e1 * t21 * (-t29 * mrSges(6,1) + t30 * mrSges(6,2)) + 0.2e1 * t39 * t99 + 0.2e1 * t48 * Ifges(5,2) * t34 + 0.2e1 * t49 * t35 * Ifges(5,1) + 0.2e1 * t47 * (-t48 * mrSges(5,1) + t49 * mrSges(5,2)) + 0.2e1 * t59 * t97 - 0.2e1 * t71 * Ifges(4,2) * t53 + 0.2e1 * t72 * t52 * Ifges(4,1) + (t53 * mrSges(4,1) + t52 * mrSges(4,2)) * t110 + (t14 * t25 + t15 * t26 + t47 * t59) * t112 + (t2 * t7 + t21 * t39 + t3 * t6) * t111 + 0.2e1 * m(4) * (t55 * t37 + t54 * t38) + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t92) * t100 + (m(4) * pkin(2) * t110 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t71 + mrSges(4,2) * t72) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t89 + (Ifges(3,1) - Ifges(3,2)) * t100) * t89) * qJD(2) + 0.2e1 * (t29 * t12 + t30 * t13) * Ifges(6,4) + 0.2e1 * (t49 * t34 + t48 * t35) * Ifges(5,4) + 0.2e1 * (t71 * t52 - t72 * t53) * Ifges(4,4) + 0.2e1 * (-t6 * t12 + t7 * t13 + t2 * t29 - t3 * t30) * mrSges(6,3) + 0.2e1 * (-t14 * t49 + t15 * t48 - t25 * t35 + t26 * t34) * mrSges(5,3) + 0.2e1 * (t37 * t71 - t38 * t72 - t52 * t54 - t53 * t55) * mrSges(4,3); t93 + (m(4) * (t37 * t88 + t38 * t91 + (-t54 * t88 + t55 * t91) * qJD(3)) + (-t91 * t52 - t88 * t53 + (t71 * t91 + t72 * t88) * qJD(3)) * mrSges(4,3)) * pkin(2) + (Ifges(3,5) * t92 - Ifges(3,6) * t89 + (-mrSges(3,1) * t92 + mrSges(3,2) * t89) * pkin(6)) * qJD(2) + m(5) * (t14 * t63 + t15 * t65 + t25 * t61 + t26 * t62) + m(6) * (t2 * t44 + t27 * t7 + t28 * t6 + t3 * t43) + (-t12 * t43 + t13 * t44 + t27 * t29 - t28 * t30) * mrSges(6,3) + (t34 * t65 - t35 * t63 + t48 * t62 - t49 * t61) * mrSges(5,3); 0.2e1 * t106 - 0.2e1 * t105 - 0.2e1 * t107 + 0.2e1 * t24 + 0.2e1 * t94 + (t61 * t63 + t62 * t65) * t112 + (t27 * t44 + t28 * t43) * t111; t93 + (m(5) * (t14 * t86 + t15 * t85) + (t34 * t85 - t35 * t86) * mrSges(5,3)) * pkin(3) + m(6) * (t2 * t66 + t3 * t64 + t57 * t7 - t58 * t6) + (-t12 * t64 + t13 * t66 + t29 * t57 + t30 * t58) * mrSges(6,3); t106 - t105 + t24 - t56 + t94 + (-t27 - t57) * mrSges(6,2) + m(6) * (t27 * t66 + t28 * t64 - t58 * t43 + t57 * t44) + m(5) * (t61 * t86 + t62 * t85) * pkin(3); 0.2e1 * m(6) * (t57 * t66 - t58 * t64) + 0.2e1 * t96; m(5) * t47 + m(6) * t21 + t97 + t99; 0; 0; 0; t95; t24 - t107; t96; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

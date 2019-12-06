% Calculate time derivative of joint inertia matrix for
% S5RRPRR5
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:38
% EndTime: 2019-12-05 18:33:40
% DurationCPUTime: 0.94s
% Computational Cost: add. (2190->182), mult. (4663->254), div. (0->0), fcn. (4289->8), ass. (0->96)
t134 = 2 * mrSges(6,3);
t84 = sin(pkin(9));
t85 = cos(pkin(9));
t109 = t84 ^ 2 + t85 ^ 2;
t108 = pkin(1) * qJD(2);
t91 = cos(qJ(2));
t99 = t91 * t108;
t77 = qJD(3) + t99;
t98 = t109 * t77;
t133 = 2 * mrSges(4,3);
t132 = 2 * mrSges(5,3);
t90 = cos(qJ(4));
t111 = t85 * t90;
t87 = sin(qJ(4));
t71 = -t84 * t87 + t111;
t72 = t84 * t90 + t85 * t87;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t44 = t71 * t89 - t72 * t86;
t45 = t71 * t86 + t72 * t89;
t24 = -mrSges(6,1) * t44 + mrSges(6,2) * t45;
t131 = 0.2e1 * t24;
t118 = pkin(8) * t72;
t88 = sin(qJ(2));
t79 = pkin(1) * t88 + qJ(3);
t62 = (-pkin(7) - t79) * t84;
t81 = t85 * pkin(7);
t63 = t79 * t85 + t81;
t41 = t90 * t62 - t63 * t87;
t29 = t41 - t118;
t42 = t87 * t62 + t90 * t63;
t65 = t71 * pkin(8);
t30 = t65 + t42;
t12 = t29 * t89 - t30 * t86;
t106 = qJD(4) * t90;
t25 = t62 * t106 + t77 * t111 + (-qJD(4) * t63 - t77 * t84) * t87;
t61 = t72 * qJD(4);
t59 = t61 * pkin(8);
t16 = -t59 + t25;
t60 = t71 * qJD(4);
t119 = pkin(8) * t60;
t26 = -t42 * qJD(4) - t72 * t77;
t17 = t26 - t119;
t2 = qJD(5) * t12 + t16 * t89 + t17 * t86;
t13 = t29 * t86 + t30 * t89;
t3 = -qJD(5) * t13 - t16 * t86 + t17 * t89;
t130 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t74 = (-pkin(7) - qJ(3)) * t84;
t76 = qJ(3) * t85 + t81;
t47 = t90 * t74 - t76 * t87;
t37 = t47 - t118;
t48 = t87 * t74 + t90 * t76;
t38 = t65 + t48;
t14 = t37 * t89 - t38 * t86;
t35 = t74 * t106 + qJD(3) * t111 + (-qJD(3) * t84 - qJD(4) * t76) * t87;
t27 = -t59 + t35;
t36 = -t72 * qJD(3) - t48 * qJD(4);
t28 = t36 - t119;
t5 = qJD(5) * t14 + t27 * t89 + t28 * t86;
t15 = t37 * t86 + t38 * t89;
t6 = -qJD(5) * t15 - t27 * t86 + t28 * t89;
t129 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t128 = -mrSges(4,1) * t85 - mrSges(5,1) * t71 + mrSges(4,2) * t84 + mrSges(5,2) * t72;
t95 = t109 * qJD(3);
t127 = 2 * m(4);
t126 = 2 * m(5);
t125 = 2 * m(6);
t22 = qJD(5) * t44 + t60 * t89 - t61 * t86;
t23 = -qJD(5) * t45 - t60 * t86 - t61 * t89;
t9 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t124 = 0.2e1 * t9;
t112 = t61 * mrSges(5,1);
t56 = t60 * mrSges(5,2);
t43 = t56 + t112;
t121 = 0.2e1 * t43;
t120 = pkin(1) * t91;
t115 = t61 * pkin(4);
t114 = t22 * mrSges(6,3);
t113 = t23 * mrSges(6,3);
t110 = Ifges(6,5) * t22 + Ifges(6,6) * t23;
t107 = pkin(4) * qJD(5);
t105 = qJD(5) * t86;
t104 = qJD(5) * t89;
t102 = t89 * t114;
t100 = t88 * t108;
t80 = -t85 * pkin(3) - pkin(2);
t97 = t109 * qJ(3);
t94 = t56 + t9;
t51 = -t71 * pkin(4) + t80;
t93 = t86 * pkin(4) * t113 + Ifges(5,5) * t60 - Ifges(5,6) * t61 + t110 + (t44 * t89 + t45 * t86) * mrSges(6,3) * t107;
t92 = 0.2e1 * t72 * t60 * Ifges(5,1) + 0.2e1 * t22 * t45 * Ifges(6,1) - 0.2e1 * t71 * Ifges(5,2) * t61 + 0.2e1 * t23 * Ifges(6,2) * t44 + 0.2e1 * (t22 * t44 + t23 * t45) * Ifges(6,4) + 0.2e1 * (t71 * t60 - t72 * t61) * Ifges(5,4);
t73 = t80 - t120;
t66 = (-mrSges(6,1) * t86 - mrSges(6,2) * t89) * t107;
t50 = t100 + t115;
t49 = t51 - t120;
t1 = [t92 + t73 * t121 + t49 * t124 + t50 * t131 - 0.2e1 * t12 * t114 + 0.2e1 * t13 * t113 - 0.2e1 * mrSges(3,2) * t99 + (t25 * t42 + t26 * t41) * t126 + (t12 * t3 + t13 * t2 + t49 * t50) * t125 - 0.2e1 * (t26 * t72 + t41 * t60) * mrSges(5,3) + (t25 * t71 - t42 * t61) * t132 + (t2 * t44 - t3 * t45) * t134 + ((-pkin(2) - t120) * t127 - (2 * mrSges(3,1)) + t73 * t126 + 0.2e1 * t128) * t100 + (t79 * t127 + t133) * t98; m(4) * (-pkin(2) * t100 + t77 * t97 + t79 * t95) + t92 + ((-t26 - t36) * t72 + (t25 + t35) * t71 - (t42 + t48) * t61 + (-t41 - t47) * t60) * mrSges(5,3) + ((-t3 - t6) * t45 + (t2 + t5) * t44 + (t13 + t15) * t23 + (-t12 - t14) * t22) * mrSges(6,3) + (t98 + t95) * mrSges(4,3) + m(5) * (t80 * t100 + t25 * t48 + t26 * t47 + t35 * t42 + t36 * t41) + m(6) * (t49 * t115 + t12 * t6 + t13 * t5 + t14 * t3 + t15 * t2 + t50 * t51) + (t49 + t51) * t9 + (t73 + t80) * t43 + (t50 + t115) * t24 + (-mrSges(3,2) * t91 + (-mrSges(3,1) + t128) * t88) * t108; (t51 * t115 + t14 * t6 + t15 * t5) * t125 + (t35 * t48 + t36 * t47) * t126 + t97 * t127 * qJD(3) + t92 + t80 * t121 + t51 * t124 + t115 * t131 + t95 * t133 + (-t14 * t22 + t15 * t23 + t5 * t44 - t6 * t45) * t134 + (t35 * t71 - t36 * t72 - t47 * t60 - t48 * t61) * t132; m(6) * t50 + t112 + (m(4) + m(5)) * t100 + t94; -(-m(6) * pkin(4) - mrSges(5,1)) * t61 + t94; 0; t26 * mrSges(5,1) - t25 * mrSges(5,2) + (-t102 + m(6) * (t13 * t104 - t12 * t105 + t2 * t86 + t3 * t89)) * pkin(4) + t93 + t130; t36 * mrSges(5,1) - t35 * mrSges(5,2) + (-t102 + m(6) * (t15 * t104 - t14 * t105 + t5 * t86 + t6 * t89)) * pkin(4) + t93 + t129; 0; 0.2e1 * t66; t110 + t130; t110 + t129; 0; t66; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

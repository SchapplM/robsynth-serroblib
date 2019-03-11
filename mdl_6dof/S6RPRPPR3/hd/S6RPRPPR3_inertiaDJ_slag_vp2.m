% Calculate time derivative of joint inertia matrix for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:08
% EndTime: 2019-03-09 02:44:10
% DurationCPUTime: 0.89s
% Computational Cost: add. (883->201), mult. (1759->278), div. (0->0), fcn. (1156->6), ass. (0->94)
t53 = sin(qJ(6));
t54 = sin(qJ(3));
t56 = cos(qJ(3));
t79 = qJD(6) * t56;
t55 = cos(qJ(6));
t84 = qJD(3) * t55;
t60 = t53 * t79 + t54 * t84;
t100 = -cos(pkin(9)) * pkin(1) - pkin(2);
t116 = 0.2e1 * t100;
t106 = -t55 / 0.2e1;
t34 = sin(pkin(9)) * pkin(1) + pkin(7);
t87 = -qJ(5) + t34;
t99 = mrSges(7,3) * t56;
t26 = -mrSges(7,2) * t54 + t53 * t99;
t27 = mrSges(7,1) * t54 + t55 * t99;
t115 = (t53 * t26 + t55 * t27) * qJD(6);
t21 = t87 * t54;
t19 = -pkin(3) * t56 - qJ(4) * t54 + t100;
t15 = pkin(4) * t56 - t19;
t9 = pkin(5) * t54 + pkin(8) * t56 + t15;
t3 = -t21 * t53 + t55 * t9;
t4 = t21 * t55 + t53 * t9;
t83 = qJD(3) * t56;
t14 = -qJD(5) * t54 + t83 * t87;
t104 = -pkin(4) - pkin(8);
t85 = qJD(3) * t54;
t42 = t54 * qJD(4);
t90 = qJ(4) * t83 + t42;
t18 = -pkin(3) * t85 + t90;
t8 = (pkin(5) * t56 + t104 * t54) * qJD(3) + t18;
t1 = qJD(6) * t3 + t14 * t55 + t53 * t8;
t2 = -qJD(6) * t4 - t14 * t53 + t55 * t8;
t71 = t1 * t55 - t2 * t53;
t114 = (t3 * t55 + t4 * t53) * qJD(6) - t71;
t113 = 2 * m(7);
t112 = 2 * mrSges(6,1);
t111 = -2 * mrSges(6,3);
t110 = 0.2e1 * t18;
t109 = -0.2e1 * t19;
t22 = t87 * t56;
t108 = 0.2e1 * t22;
t107 = t53 / 0.2e1;
t105 = m(5) + m(6);
t103 = m(5) * t18;
t52 = qJ(4) + pkin(5);
t102 = m(7) * t52;
t101 = t3 * t53;
t98 = Ifges(7,4) * t53;
t97 = Ifges(7,4) * t55;
t96 = Ifges(7,5) * t55;
t13 = qJD(5) * t56 + t85 * t87;
t95 = t13 * t22;
t65 = Ifges(7,2) * t55 + t98;
t94 = t53 * t65;
t67 = Ifges(7,1) * t53 + t97;
t93 = t55 * t67;
t92 = -mrSges(4,1) - mrSges(5,1);
t91 = mrSges(4,2) - mrSges(5,3);
t89 = mrSges(6,1) * t83 + mrSges(6,2) * t85;
t88 = t53 ^ 2 + t55 ^ 2;
t86 = qJD(3) * t22;
t82 = qJD(4) * t22;
t81 = qJD(6) * t53;
t80 = qJD(6) * t55;
t77 = t55 * t79;
t76 = t53 * t85;
t74 = m(7) * t88;
t73 = Ifges(7,5) * t60 + Ifges(7,6) * t77 + Ifges(7,3) * t83;
t72 = m(5) * t34 + mrSges(5,2) - mrSges(6,3);
t28 = mrSges(7,1) * t55 - mrSges(7,2) * t53;
t69 = -mrSges(7,1) * t53 - mrSges(7,2) * t55;
t68 = Ifges(7,1) * t55 - t98;
t66 = -Ifges(7,2) * t53 + t97;
t10 = mrSges(7,1) * t83 - mrSges(7,3) * t60;
t59 = t76 - t77;
t11 = -mrSges(7,2) * t83 - mrSges(7,3) * t59;
t64 = t53 * t10 - t55 * t11;
t63 = t26 * t55 - t27 * t53;
t61 = -Ifges(7,6) * t53 - (2 * Ifges(4,4)) - (2 * Ifges(6,4)) + (2 * Ifges(5,5));
t58 = -m(7) * t114 - t64;
t57 = -pkin(3) - pkin(4);
t48 = -pkin(3) + t104;
t40 = Ifges(7,6) * t81;
t25 = t68 * qJD(6);
t24 = t66 * qJD(6);
t23 = t69 * qJD(6);
t20 = t69 * t56;
t17 = Ifges(7,5) * t54 - t56 * t68;
t16 = Ifges(7,6) * t54 - t56 * t66;
t12 = -pkin(4) * t85 + t18;
t7 = mrSges(7,1) * t59 + mrSges(7,2) * t60;
t6 = t67 * t79 + (t56 * Ifges(7,5) + t54 * t68) * qJD(3);
t5 = t65 * t79 + (t56 * Ifges(7,6) + t54 * t66) * qJD(3);
t29 = [0.2e1 * t15 * t89 + t103 * t109 + 0.2e1 * t1 * t26 + 0.2e1 * t2 * t27 - 0.2e1 * t13 * t20 + t7 * t108 + 0.2e1 * t3 * t10 + 0.2e1 * t4 * t11 + (t1 * t4 + t2 * t3 - t95) * t113 + 0.2e1 * m(6) * (t12 * t15 + t14 * t21 - t95) + (mrSges(5,3) * t110 + t14 * t111 + t12 * t112 + t73) * t54 + (mrSges(5,1) * t110 - 0.2e1 * t12 * mrSges(6,2) + 0.2e1 * t13 * mrSges(6,3) + t53 * t5 - t55 * t6 + (t16 * t55 + t17 * t53) * qJD(6)) * t56 + ((t21 * t111 + mrSges(5,3) * t109 + mrSges(4,2) * t116 + (-t61 - t96) * t56) * t56 + (0.2e1 * t19 * mrSges(5,1) + mrSges(4,1) * t116 - t53 * t16 + mrSges(6,3) * t108 + t55 * t17 + t61 * t54 + ((2 * Ifges(6,2)) + (2 * Ifges(5,1)) - (2 * Ifges(5,3)) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(6,1)) + Ifges(7,3)) * t56) * t54) * qJD(3); (t7 + t63 * qJD(3) + m(7) * (-qJD(3) * t101 + t4 * t84 - t13) + m(6) * (qJD(3) * t21 - t13)) * t54 + (qJD(3) * t20 + t115 + m(7) * (t3 * t80 + t4 * t81 - t71 + t86) + m(6) * (-t14 + t86) + t64) * t56; (0.1e1 - t88) * t54 * t83 * t113; t54 * (-Ifges(7,5) * t80 + t40) / 0.2e1 + t5 * t106 + t52 * t7 - t17 * t80 / 0.2e1 + t22 * t23 + t14 * mrSges(6,2) + qJD(4) * t20 + (-t6 / 0.2e1 + qJD(6) * t16 / 0.2e1) * t53 + (-t28 - mrSges(6,1)) * t13 + m(6) * (-qJ(4) * t13 + t14 * t57 + t82) + m(7) * (-t13 * t52 + t82) + t114 * mrSges(7,3) + (-t26 * t81 - t27 * t80 + t58) * t48 + (-t25 * t106 - t24 * t107 + (t106 * t65 - t107 * t67) * qJD(6) + t72 * qJD(4)) * t56 + ((-pkin(3) * mrSges(5,2) - t57 * mrSges(6,3) - Ifges(7,5) * t53 / 0.2e1 + Ifges(7,6) * t106 + Ifges(6,6) + Ifges(5,4) + Ifges(4,5) + (-m(5) * pkin(3) + t92) * t34) * t56 + (-Ifges(6,5) + Ifges(5,6) - Ifges(4,6) + t94 / 0.2e1 - t93 / 0.2e1 + t91 * t34 - t72 * qJ(4)) * t54) * qJD(3); t54 * t23 + m(7) * t42 + t103 + m(6) * t90 + ((t28 - t91 + t102) * t56 + (m(6) * t57 - mrSges(7,3) * t88 + t48 * t74 + t92) * t54) * qJD(3) + t89; 0.2e1 * t23 * t52 + t24 * t55 + t25 * t53 + (t93 - t94) * qJD(6) + (0.2e1 * qJ(4) * t105 + 0.2e1 * mrSges(5,3) + 0.2e1 * t102 + t112 + 0.2e1 * t28) * qJD(4); m(6) * t14 + t72 * t83 - t115 + t58; (t74 + t105) * t85; 0; 0; t55 * t10 + t53 * t11 + t63 * qJD(6) + m(7) * (t1 * t53 + t2 * t55 + (t4 * t55 - t101) * qJD(6)) + m(6) * t12 + t89; 0; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t76 + t73; -t7; t40 + (-t28 * t48 - t96) * qJD(6); -t28 * qJD(6); t23; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t29(1) t29(2) t29(4) t29(7) t29(11) t29(16); t29(2) t29(3) t29(5) t29(8) t29(12) t29(17); t29(4) t29(5) t29(6) t29(9) t29(13) t29(18); t29(7) t29(8) t29(9) t29(10) t29(14) t29(19); t29(11) t29(12) t29(13) t29(14) t29(15) t29(20); t29(16) t29(17) t29(18) t29(19) t29(20) t29(21);];
Mq  = res;

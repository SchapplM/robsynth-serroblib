% Calculate time derivative of joint inertia matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:14
% EndTime: 2019-03-09 01:37:16
% DurationCPUTime: 1.21s
% Computational Cost: add. (1044->223), mult. (2052->356), div. (0->0), fcn. (1551->6), ass. (0->106)
t53 = sin(qJ(6));
t55 = cos(qJ(6));
t36 = -mrSges(7,1) * t55 + mrSges(7,2) * t53;
t70 = -m(7) * pkin(5) - mrSges(6,1) + t36;
t49 = sin(pkin(9));
t50 = cos(pkin(9));
t35 = qJD(2) * t49 - qJD(3) * t50;
t117 = 0.2e1 * t35;
t89 = t53 ^ 2 + t55 ^ 2;
t116 = m(7) * pkin(8) + mrSges(7,3);
t66 = mrSges(7,1) * t53 + mrSges(7,2) * t55;
t28 = t66 * qJD(6);
t54 = sin(qJ(5));
t56 = cos(qJ(5));
t86 = qJD(5) * t56;
t115 = t54 * (qJD(5) * mrSges(6,2) + t28) + t70 * t86;
t114 = 0.2e1 * m(7);
t51 = qJ(2) + pkin(3);
t52 = -pkin(1) - qJ(3);
t90 = t49 * t51 + t50 * t52;
t27 = pkin(7) + t90;
t113 = 0.2e1 * t27;
t34 = t50 * qJD(2) + t49 * qJD(3);
t112 = 0.2e1 * t34;
t106 = Ifges(7,4) * t53;
t37 = Ifges(7,2) * t55 + t106;
t111 = -t37 / 0.2e1;
t110 = -t53 / 0.2e1;
t109 = pkin(5) * t54;
t108 = pkin(8) * t56;
t107 = mrSges(7,3) * t54;
t105 = Ifges(7,4) * t55;
t104 = Ifges(7,5) * t53;
t103 = Ifges(7,6) * t53;
t102 = Ifges(7,6) * t55;
t101 = Ifges(7,6) * t56;
t48 = t56 ^ 2;
t100 = t27 * t48;
t99 = t34 * t49;
t98 = t34 * t50;
t97 = t35 * t49;
t96 = t35 * t50;
t95 = t35 * t56;
t46 = t54 ^ 2;
t94 = t46 * t35;
t93 = t53 * t56;
t92 = t55 * t56;
t75 = t55 * t86;
t87 = qJD(5) * t54;
t91 = -Ifges(7,5) * t75 - Ifges(7,3) * t87;
t88 = t46 - t48;
t85 = qJD(6) * t53;
t84 = qJD(6) * t54;
t83 = qJD(6) * t55;
t81 = t53 * t87;
t80 = t55 * t87;
t79 = t54 * t86;
t76 = t53 * t84;
t74 = (2 * Ifges(6,4)) + t103;
t72 = -t49 * t52 + t50 * t51;
t26 = -pkin(4) - t72;
t15 = -pkin(5) * t56 - t54 * pkin(8) + t26;
t59 = -qJD(6) * t15 + t27 * t87 - t95;
t68 = -qJD(6) * t27 * t56 + (-t108 + t109) * qJD(5) - t34;
t1 = t68 * t53 - t59 * t55;
t3 = t55 * t15 - t27 * t93;
t73 = -qJD(6) * t3 + t1;
t71 = 0.2e1 * t79;
t69 = t50 * t71;
t67 = t54 * mrSges(6,1) + t56 * mrSges(6,2);
t65 = Ifges(7,1) * t55 - t106;
t38 = Ifges(7,1) * t53 + t105;
t64 = -Ifges(7,2) * t53 + t105;
t24 = t53 * t49 + t50 * t92;
t22 = t55 * t49 - t50 * t93;
t23 = t49 * t92 - t53 * t50;
t21 = -t49 * t93 - t55 * t50;
t25 = t66 * t54;
t60 = t53 * t86 + t54 * t83;
t61 = t75 - t76;
t7 = t60 * mrSges(7,1) + t61 * mrSges(7,2);
t62 = t25 * t86 + t54 * t7;
t10 = -t23 * qJD(6) + t49 * t81;
t8 = t21 * qJD(6) - t49 * t80;
t58 = -t10 * t53 + t8 * t55 + (-t21 * t55 - t23 * t53) * qJD(6);
t11 = -t24 * qJD(6) + t50 * t81;
t9 = t22 * qJD(6) - t50 * t80;
t57 = -t11 * t53 + t9 * t55 + (-t22 * t55 - t24 * t53) * qJD(6);
t44 = Ifges(7,5) * t83;
t33 = -mrSges(7,1) * t56 - t55 * t107;
t32 = mrSges(7,2) * t56 - t53 * t107;
t31 = t65 * qJD(6);
t30 = t64 * qJD(6);
t29 = t67 * qJD(5);
t20 = t50 * t94;
t19 = t49 * t94;
t18 = -t56 * Ifges(7,5) + t65 * t54;
t17 = t64 * t54 - t101;
t14 = -mrSges(7,2) * t87 - t60 * mrSges(7,3);
t13 = mrSges(7,1) * t87 - t61 * mrSges(7,3);
t12 = t27 * t94;
t6 = -t38 * t84 + (Ifges(7,5) * t54 + t65 * t56) * qJD(5);
t5 = -t37 * t84 + (Ifges(7,6) * t54 + t64 * t56) * qJD(5);
t4 = t53 * t15 + t27 * t92;
t2 = t59 * t53 + t68 * t55;
t16 = [mrSges(5,1) * t112 + 0.2e1 * qJD(3) * mrSges(4,3) + 0.2e1 * t1 * t32 + 0.2e1 * t3 * t13 + 0.2e1 * t4 * t14 + 0.2e1 * t2 * t33 + 0.2e1 * t26 * t29 - 0.2e1 * t35 * mrSges(5,2) + (t4 * t1 + t3 * t2 + t12) * t114 + 0.2e1 * m(5) * (t72 * t34 + t90 * t35) + 0.2e1 * m(4) * (qJ(2) * qJD(2) - qJD(3) * t52) + 0.2e1 * m(6) * (t35 * t100 - t26 * t34 + t12) + (mrSges(6,1) * t112 + (t25 * t113 - t53 * t17 + t55 * t18 + t74 * t56) * qJD(5) + t91) * t56 + (-0.2e1 * t34 * mrSges(6,2) + t25 * t117 + t7 * t113 - t53 * t5 + t55 * t6 + (-t55 * t17 - t53 * t18 - t56 * (-t102 - t104)) * qJD(6) + ((Ifges(7,5) * t55 - t74) * t54 + (t27 ^ 2 * t114 + (2 * Ifges(6,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t56) * qJD(5)) * t54 + (t46 + t48) * mrSges(6,3) * t117 + 0.2e1 * (m(3) * qJ(2) + mrSges(4,1) + mrSges(3,3)) * qJD(2); -m(4) * qJD(3) + t11 * t33 + t22 * t13 + t24 * t14 + t49 * t29 + t9 * t32 + t62 * t50 + m(7) * (t24 * t1 + t11 * t3 + t22 * t2 + t27 * t69 + t9 * t4 + t20) + m(6) * (t48 * t96 + t20 - t99) + m(5) * (t96 - t99); (t50 ^ 2 * t79 + t22 * t11 + t24 * t9) * t114; m(4) * qJD(2) + t10 * t33 + t21 * t13 + t23 * t14 - t50 * t29 + t8 * t32 + t62 * t49 + m(7) * (t27 * t49 * t71 + t23 * t1 + t10 * t3 + t21 * t2 + t8 * t4 + t19) + m(6) * (t48 * t97 + t19 + t98) + m(5) * (t97 + t98); m(7) * (t10 * t22 + t21 * t11 + t23 * t9 + t8 * t24 + t49 * t69); (t49 ^ 2 * t79 + t21 * t10 + t23 * t8) * t114; -t56 * t7 + (t32 * t92 - t33 * t93 + m(7) * (t27 * t46 - t3 * t93 + t4 * t92 - t100)) * qJD(5) + (-t32 * t85 + t55 * t14 - t33 * t83 - t53 * t13 + m(7) * (t1 * t55 - t2 * t53 - t3 * t83 - t4 * t85 - t95) + qJD(5) * t25) * t54; m(7) * (t57 * t54 + ((-t22 * t53 + t24 * t55) * t56 + t88 * t50) * qJD(5)); m(7) * (t58 * t54 + ((-t21 * t53 + t23 * t55) * t56 + t88 * t49) * qJD(5)); m(7) * (-0.1e1 + t89) * t71; -pkin(5) * t7 + (-t35 * mrSges(6,2) - t44 / 0.2e1 + (t70 * t27 + Ifges(6,5)) * qJD(5)) * t56 + (t86 * t111 - t2 * mrSges(7,3) + t6 / 0.2e1 + (-t4 * mrSges(7,3) - t17 / 0.2e1 + t101 / 0.2e1) * qJD(6) + (-qJD(6) * t32 - t13 + m(7) * (-qJD(6) * t4 - t2)) * pkin(8)) * t53 + (t38 * t86 / 0.2e1 + qJD(6) * t18 / 0.2e1 + t5 / 0.2e1 + t73 * mrSges(7,3) + (m(7) * t73 - qJD(6) * t33 + t14) * pkin(8)) * t55 + (t55 * t31 / 0.2e1 + t30 * t110 + t27 * t28 + (t38 * t110 + t55 * t111) * qJD(6) + (t27 * mrSges(6,2) + t104 / 0.2e1 + t102 / 0.2e1 - Ifges(6,6)) * qJD(5) + t70 * t35) * t54; t115 * t50 + t116 * t57; t115 * t49 + t116 * t58; -t56 * t28 + (m(7) * (t89 * t108 - t109) + t54 * t36 + t89 * t56 * mrSges(7,3) - t67) * qJD(5); -0.2e1 * pkin(5) * t28 + t30 * t55 + t31 * t53 + (-t37 * t53 + t38 * t55) * qJD(6); t2 * mrSges(7,1) - t1 * mrSges(7,2) - Ifges(7,5) * t76 - t60 * Ifges(7,6) - t91; mrSges(7,1) * t11 - mrSges(7,2) * t9; mrSges(7,1) * t10 - mrSges(7,2) * t8; -t7; t44 + (t36 * pkin(8) - t103) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;

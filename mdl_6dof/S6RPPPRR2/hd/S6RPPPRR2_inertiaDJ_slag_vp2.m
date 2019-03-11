% Calculate time derivative of joint inertia matrix for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:40
% EndTime: 2019-03-09 01:31:42
% DurationCPUTime: 0.86s
% Computational Cost: add. (1528->203), mult. (3004->309), div. (0->0), fcn. (2702->8), ass. (0->90)
t107 = 2 * qJD(3);
t53 = sin(pkin(10));
t55 = cos(pkin(10));
t68 = (t53 ^ 2 + t55 ^ 2) * qJD(4);
t96 = sin(qJ(5));
t71 = qJD(5) * t96;
t97 = cos(qJ(5));
t72 = qJD(5) * t97;
t29 = -t53 * t71 + t55 * t72;
t31 = t97 * t53 + t96 * t55;
t21 = t31 * t29;
t28 = -t53 * t72 - t55 * t71;
t30 = t96 * t53 - t97 * t55;
t22 = t30 * t28;
t106 = m(6) * (t21 - t22);
t56 = sin(qJ(6));
t57 = cos(qJ(6));
t36 = -mrSges(7,1) * t57 + t56 * mrSges(7,2);
t105 = m(7) * pkin(5) + mrSges(6,1) - t36;
t104 = -2 * mrSges(6,3);
t43 = -cos(pkin(9)) * pkin(1) - pkin(2) - qJ(4);
t98 = -pkin(7) + t43;
t27 = t98 * t53;
t67 = t98 * t97;
t15 = t96 * t27 - t55 * t67;
t103 = 0.2e1 * t15;
t102 = t30 / 0.2e1;
t93 = Ifges(7,4) * t57;
t38 = Ifges(7,1) * t56 + t93;
t101 = t38 / 0.2e1;
t100 = pkin(5) * t29;
t66 = t98 * t96;
t69 = t96 * qJD(4);
t70 = t97 * qJD(4);
t9 = t27 * t72 - t53 * t69 + (qJD(5) * t66 + t70) * t55;
t99 = t15 * t9;
t95 = mrSges(7,3) * t30;
t94 = Ifges(7,4) * t56;
t92 = Ifges(7,5) * t29;
t91 = Ifges(7,6) * t29;
t90 = Ifges(7,6) * t31;
t89 = Ifges(7,6) * t56;
t88 = t57 * t28;
t87 = Ifges(7,5) * t88 + Ifges(7,3) * t29;
t85 = t56 ^ 2 + t57 ^ 2;
t44 = sin(pkin(9)) * pkin(1) + qJ(3);
t35 = t53 * pkin(4) + t44;
t14 = pkin(5) * t31 + pkin(8) * t30 + t35;
t16 = t97 * t27 + t55 * t66;
t5 = t14 * t57 - t56 * t16;
t84 = t5 * qJD(6);
t6 = t56 * t14 + t16 * t57;
t83 = t6 * qJD(6);
t82 = qJD(6) * t30;
t81 = qJD(6) * t56;
t80 = qJD(6) * t57;
t79 = t29 * t104;
t77 = t30 * t81;
t76 = t85 * t28;
t75 = t85 * t29;
t74 = t29 * mrSges(6,1) + t28 * mrSges(6,2);
t73 = -(2 * Ifges(6,4)) - t89;
t65 = -t5 * t56 + t57 * t6;
t64 = -t28 * t15 + t30 * t9;
t63 = t15 * t29 + t9 * t31;
t62 = mrSges(7,1) * t56 + mrSges(7,2) * t57;
t61 = Ifges(7,1) * t57 - t94;
t60 = -Ifges(7,2) * t56 + t93;
t59 = -t56 * t28 + t30 * t80;
t58 = t77 + t88;
t46 = Ifges(7,5) * t80;
t37 = Ifges(7,2) * t57 + t94;
t34 = t61 * qJD(6);
t33 = t60 * qJD(6);
t32 = t62 * qJD(6);
t20 = t31 * mrSges(7,1) + t57 * t95;
t19 = -mrSges(7,2) * t31 + t56 * t95;
t18 = t62 * t30;
t17 = -pkin(8) * t28 + qJD(3) + t100;
t13 = Ifges(7,5) * t31 - t61 * t30;
t12 = -t60 * t30 + t90;
t11 = -t29 * mrSges(7,2) + t59 * mrSges(7,3);
t10 = t29 * mrSges(7,1) - t58 * mrSges(7,3);
t8 = -t53 * t70 - t27 * t71 + (qJD(5) * t67 - t69) * t55;
t7 = -t59 * mrSges(7,1) + t58 * mrSges(7,2);
t4 = t58 * Ifges(7,1) + t59 * Ifges(7,4) + t92;
t3 = t58 * Ifges(7,4) + t59 * Ifges(7,2) + t91;
t2 = t17 * t57 - t56 * t8 - t83;
t1 = t56 * t17 + t57 * t8 + t84;
t23 = [t16 * t79 + 0.2e1 * t35 * t74 + t7 * t103 - 0.2e1 * t9 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t5 * t10 + 0.2e1 * t6 * t11 + (mrSges(6,3) * t103 - t56 * t12 + t57 * t13) * t28 + 0.2e1 * m(5) * (qJD(3) * t44 - t43 * t68) + 0.2e1 * m(6) * (qJD(3) * t35 + t16 * t8 + t99) + 0.2e1 * m(7) * (t1 * t6 + t2 * t5 + t99) + (mrSges(6,1) * t107 + t8 * t104 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t29 + t73 * t28 + t87) * t31 + (-0.2e1 * qJD(3) * mrSges(6,2) + t9 * t104 - 0.2e1 * Ifges(6,1) * t28 + t56 * t3 - t57 * t4 + (-Ifges(7,5) * t57 - t73) * t29 + (t31 * (Ifges(7,5) * t56 + Ifges(7,6) * t57) + t56 * t13 + t57 * t12) * qJD(6)) * t30 + (m(4) * t44 + mrSges(5,1) * t53 + mrSges(5,2) * t55 + mrSges(4,3)) * t107 + 0.2e1 * mrSges(5,3) * t68; -t29 * t18 + t31 * t7 + m(6) * (t16 * t28 - t30 * t8 + t63) + m(7) * t63 + (t20 * t82 + t28 * t19 - t30 * t11 + m(7) * (-t1 * t30 + t28 * t6 + t5 * t82)) * t57 + (-t28 * t20 + t30 * t10 + t19 * t82 + m(7) * (t2 * t30 - t28 * t5 + t6 * t82)) * t56; 0.2e1 * t106 + 0.2e1 * m(7) * (-t85 * t22 + t21); t30 * t7 + (t57 * t19 - t56 * t20) * t29 + (0.2e1 * t30 * mrSges(6,3) + t18) * t28 + (t79 - t56 * t10 + t57 * t11 + (-t19 * t56 - t20 * t57) * qJD(6)) * t31 + m(7) * (t65 * t29 + (t1 * t57 - t2 * t56 + (-t5 * t57 - t56 * t6) * qJD(6)) * t31 + t64) + m(6) * (t16 * t29 + t8 * t31 + t64) - m(5) * t68; m(7) * (t28 * t31 - t30 * t29) * (-0.1e1 + t85); 0.2e1 * t106 + 0.2e1 * m(7) * (t31 * t75 - t22); t19 * t80 + t56 * t11 - t20 * t81 + t57 * t10 + m(7) * (t65 * qJD(6) + t56 * t1 + t2 * t57) + (m(6) + m(5)) * qJD(3) + t74; 0; 0; 0; t31 * t46 / 0.2e1 + t15 * t32 - Ifges(6,6) * t29 + Ifges(6,5) * t28 - pkin(5) * t7 - t8 * mrSges(6,2) - t105 * t9 + (t33 * t102 - t28 * t37 / 0.2e1 - t2 * mrSges(7,3) + t4 / 0.2e1 + t92 / 0.2e1 + (-t12 / 0.2e1 - t90 / 0.2e1 + t30 * t101 - t6 * mrSges(7,3)) * qJD(6) + (-t10 + m(7) * (-t2 - t83) - qJD(6) * t19) * pkin(8)) * t56 + (t28 * t101 + t1 * mrSges(7,3) - t30 * t34 / 0.2e1 + t91 / 0.2e1 + t3 / 0.2e1 + (t13 / 0.2e1 + t37 * t102 - t5 * mrSges(7,3)) * qJD(6) + (t11 - qJD(6) * t20 + m(7) * (t1 - t84)) * pkin(8)) * t57; t29 * t36 + t31 * t32 + m(7) * (pkin(8) * t76 - t100) + mrSges(7,3) * t76 - t74; -t29 * mrSges(6,2) + t30 * t32 + (m(7) * pkin(8) + mrSges(7,3)) * t75 + t105 * t28; 0; -0.2e1 * pkin(5) * t32 + t33 * t57 + t56 * t34 + (-t56 * t37 + t57 * t38) * qJD(6); t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t77 + t59 * Ifges(7,6) + t87; -t7; (-t57 * t29 + t31 * t81) * mrSges(7,2) + (-t29 * t56 - t31 * t80) * mrSges(7,1); -t32; t46 + (t36 * pkin(8) - t89) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;

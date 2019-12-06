% Calculate time derivative of joint inertia matrix for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:00
% EndTime: 2019-12-05 16:12:04
% DurationCPUTime: 0.98s
% Computational Cost: add. (484->199), mult. (1438->310), div. (0->0), fcn. (996->6), ass. (0->99)
t59 = cos(pkin(8));
t55 = t59 ^ 2;
t58 = sin(pkin(8));
t120 = qJD(4) * (t58 ^ 2 + t55);
t119 = -(Ifges(5,6) - Ifges(6,6)) * t58 + (Ifges(6,4) + Ifges(5,5)) * t59;
t118 = m(5) + m(6);
t117 = mrSges(6,2) + mrSges(5,3);
t116 = Ifges(6,1) + Ifges(5,1);
t113 = 2 * m(5);
t112 = 2 * m(6);
t111 = 2 * pkin(6);
t65 = pkin(4) * t58 - qJ(5) * t59 + pkin(6);
t60 = sin(qJ(3));
t79 = qJD(5) * t60;
t62 = cos(qJ(3));
t81 = qJD(3) * t62;
t5 = -t59 * t79 + t65 * t81;
t110 = m(6) * t5;
t61 = sin(qJ(2));
t82 = qJD(3) * t61;
t72 = t60 * t82;
t63 = cos(qJ(2));
t94 = t59 * t63;
t9 = -t59 * t72 + (t58 * t61 + t62 * t94) * qJD(2);
t108 = t59 * t9;
t92 = t61 * t62;
t25 = -t58 * t63 + t59 * t92;
t107 = t59 * qJD(4) * t25 + qJ(4) * t108;
t106 = mrSges(4,2) * t62;
t105 = mrSges(6,3) * t59;
t104 = Ifges(5,4) * t58;
t103 = Ifges(5,4) * t59;
t102 = Ifges(6,5) * t58;
t101 = Ifges(6,5) * t59;
t23 = -qJD(4) * t60 + (pkin(3) * t60 - qJ(4) * t62) * qJD(3);
t100 = t23 * t59;
t39 = -pkin(3) * t62 - qJ(4) * t60 - pkin(2);
t99 = t39 * t59;
t98 = t58 * t60;
t97 = t58 * t62;
t96 = t59 * t60;
t95 = t59 * t62;
t93 = t60 * mrSges(4,2);
t26 = (mrSges(6,1) * t58 - t105) * t60;
t27 = (mrSges(5,1) * t58 + mrSges(5,2) * t59) * t60;
t91 = t26 + t27;
t28 = (-mrSges(5,2) * t60 - mrSges(5,3) * t97) * qJD(3);
t31 = (-mrSges(6,2) * t97 + mrSges(6,3) * t60) * qJD(3);
t90 = t28 + t31;
t29 = (mrSges(5,1) * t60 - mrSges(5,3) * t95) * qJD(3);
t73 = t59 * t81;
t83 = qJD(3) * t60;
t30 = -mrSges(6,1) * t83 + mrSges(6,2) * t73;
t89 = -t29 + t30;
t33 = mrSges(5,2) * t62 - mrSges(5,3) * t98;
t36 = -mrSges(6,2) * t98 - mrSges(6,3) * t62;
t88 = t33 + t36;
t34 = -mrSges(5,1) * t62 - mrSges(5,3) * t96;
t35 = mrSges(6,1) * t62 + mrSges(6,2) * t96;
t87 = -t34 + t35;
t86 = -mrSges(5,1) * t59 + mrSges(5,2) * t58 - mrSges(4,1);
t74 = t58 * t81;
t22 = mrSges(5,1) * t74 + mrSges(5,2) * t73;
t17 = pkin(6) * t95 + t58 * t39;
t85 = qJ(4) * t120;
t84 = qJD(2) * t63;
t80 = qJD(5) * t58;
t41 = -mrSges(6,1) * t59 - mrSges(6,3) * t58;
t77 = t41 + t86;
t76 = pkin(6) * t83;
t75 = t61 * t84;
t71 = t61 * t81;
t70 = pkin(6) * t58 + pkin(4);
t24 = t58 * t92 + t94;
t8 = -qJD(2) * t61 * t59 - t58 * t72 + t84 * t97;
t66 = qJ(4) * t8 + qJD(4) * t24;
t64 = t60 * t84 + t71;
t57 = t62 ^ 2;
t56 = t60 ^ 2;
t49 = t56 * pkin(6) * t84;
t45 = mrSges(6,1) * t74;
t43 = t56 * t75;
t38 = -pkin(4) * t59 - qJ(5) * t58 - pkin(3);
t37 = (mrSges(4,1) * t60 + t106) * qJD(3);
t21 = -mrSges(6,3) * t73 + t45;
t20 = t65 * t60;
t19 = t58 * t23;
t16 = -pkin(6) * t97 + t99;
t15 = (t60 * Ifges(5,5) + (t59 * Ifges(5,1) - t104) * t62) * qJD(3);
t14 = (t60 * Ifges(6,4) + (t59 * Ifges(6,1) + t102) * t62) * qJD(3);
t13 = (t60 * Ifges(5,6) + (-t58 * Ifges(5,2) + t103) * t62) * qJD(3);
t12 = (t60 * Ifges(6,6) + (t58 * Ifges(6,3) + t101) * t62) * qJD(3);
t11 = t62 * t70 - t99;
t10 = -qJ(5) * t62 + t17;
t7 = -t59 * t76 + t19;
t6 = t58 * t76 + t100;
t3 = -t70 * t83 - t100;
t2 = -qJD(5) * t62 + t19 + (-pkin(6) * t59 + qJ(5)) * t83;
t1 = [0.2e1 * m(4) * (t43 + (t57 - 0.1e1) * t75) + 0.2e1 * t118 * (t61 ^ 2 * t60 * t81 + t24 * t8 + t25 * t9 + t43); -t63 * t37 + t88 * t9 + t87 * t8 + t90 * t25 + t89 * t24 + ((t21 + t22) * t60 + t91 * t81) * t61 + ((-t62 * mrSges(4,1) - mrSges(3,1) + t93) * t61 + (-mrSges(3,2) + t91 * t60 + (t56 + t57) * mrSges(4,3)) * t63) * qJD(2) + m(4) * (t49 + (pkin(6) * t57 * t63 - pkin(2) * t61) * qJD(2)) + m(5) * (t111 * t60 * t71 - t16 * t8 + t17 * t9 - t24 * t6 + t25 * t7 + t49) + m(6) * (t5 * t60 * t61 + t10 * t9 + t11 * t8 + t2 * t25 + t20 * t64 + t3 * t24); -0.2e1 * pkin(2) * t37 + 0.2e1 * t10 * t31 + 0.2e1 * t11 * t30 + 0.2e1 * t16 * t29 + 0.2e1 * t17 * t28 + 0.2e1 * t2 * t36 + 0.2e1 * t20 * t21 + 0.2e1 * t5 * t26 + 0.2e1 * t3 * t35 + 0.2e1 * t7 * t33 + 0.2e1 * t6 * t34 + (t16 * t6 + t17 * t7) * t113 + (t10 * t2 + t11 * t3 + t20 * t5) * t112 + (t22 * t111 + (t14 + t15) * t59 + (t12 - t13) * t58 + (-(2 * Ifges(4,4)) + t119) * t83) * t60 + (t27 * t111 + 0.2e1 * (Ifges(4,4) - t119) * t62 + (-(2 * Ifges(6,2)) - (2 * Ifges(5,3)) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + (pkin(6) ^ 2 * t113) + t116 * t55 + ((Ifges(6,3) + Ifges(5,2)) * t58 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t59) * t58) * t60) * t81; (t62 * t77 + t93) * t82 + (t60 * t77 - t106) * t84 + m(5) * (-pkin(3) * t64 + t58 * t66 + t107) + m(6) * (t64 * t38 + (-t61 * t79 + t66) * t58 + t107) + t117 * (t58 * t8 + t108); -pkin(3) * t22 + t5 * t41 + (t21 + t110) * t38 + (-t12 / 0.2e1 + t13 / 0.2e1 + t2 * mrSges(6,2) + t7 * mrSges(5,3) + t88 * qJD(4) + t90 * qJ(4) + m(5) * (qJ(4) * t7 + qJD(4) * t17) + m(6) * (qJ(4) * t2 + qJD(4) * t10)) * t59 + (t14 / 0.2e1 + t15 / 0.2e1 + t3 * mrSges(6,2) - t6 * mrSges(5,3) - qJD(5) * t26 + t87 * qJD(4) + t89 * qJ(4) + m(5) * (-qJ(4) * t6 - qJD(4) * t16) + m(6) * (qJ(4) * t3 + qJD(4) * t11 - qJD(5) * t20)) * t58 + ((pkin(6) * mrSges(4,2) - Ifges(4,6) + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t59 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t58) * t60 + (t58 * (-Ifges(6,3) * t59 + t102) / 0.2e1 - t58 * (Ifges(5,2) * t59 + t104) / 0.2e1 + Ifges(4,5) + (-m(5) * pkin(3) + t86) * pkin(6) + (t116 * t58 - t101 + t103) * t59 / 0.2e1) * t62) * qJD(3); -0.2e1 * t41 * t80 + t85 * t113 + (-t38 * t80 + t85) * t112 + 0.2e1 * t117 * t120; t118 * t64; t110 + t45 + ((m(5) * pkin(6)) - t105) * t81 + t22; -m(6) * t80; 0; m(6) * t8; m(6) * t3 + t30; m(6) * t58 * qJD(4); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

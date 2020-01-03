% Calculate time derivative of joint inertia matrix for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR15_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:39
% EndTime: 2019-12-31 18:36:42
% DurationCPUTime: 1.09s
% Computational Cost: add. (1120->228), mult. (2590->370), div. (0->0), fcn. (2090->6), ass. (0->99)
t81 = (-pkin(1) - pkin(6));
t114 = 2 * t81;
t76 = cos(pkin(8));
t74 = t76 ^ 2;
t75 = sin(pkin(8));
t95 = t75 ^ 2 + t74;
t113 = qJD(4) * m(5) * t95;
t77 = sin(qJ(5));
t79 = cos(qJ(5));
t83 = t75 * t77 - t76 * t79;
t112 = qJD(5) * t83;
t111 = 2 * m(6);
t110 = 2 * mrSges(4,1);
t109 = -t83 / 0.2e1;
t62 = t75 * t79 + t76 * t77;
t108 = t62 / 0.2e1;
t107 = t75 / 0.2e1;
t78 = sin(qJ(3));
t106 = pkin(3) * t78;
t80 = cos(qJ(3));
t105 = pkin(7) * t80;
t104 = Ifges(5,4) * t75;
t103 = Ifges(5,4) * t76;
t102 = t75 * Ifges(5,2);
t101 = t76 * t78;
t100 = t78 * t81;
t99 = t80 * mrSges(5,3);
t98 = pkin(7) + qJ(4);
t54 = t62 * qJD(5);
t97 = -Ifges(6,5) * t112 - Ifges(6,6) * t54;
t96 = -mrSges(5,1) * t76 + mrSges(5,2) * t75 - mrSges(4,1);
t48 = -qJD(4) * t80 + qJD(2) + (pkin(3) * t80 + qJ(4) * t78) * qJD(3);
t92 = qJD(3) * t80;
t90 = t81 * t92;
t29 = t75 * t48 + t76 * t90;
t94 = qJ(4) * t80;
t65 = qJ(2) - t94 + t106;
t41 = t76 * t100 + t75 * t65;
t93 = qJD(3) * t78;
t19 = -t80 * t54 + t83 * t93;
t21 = t112 * t80 + t62 * t93;
t91 = Ifges(6,5) * t19 + Ifges(6,6) * t21 + Ifges(6,3) * t92;
t88 = -t75 * t81 + pkin(4);
t87 = pkin(4) * t75 - t81;
t86 = t95 * mrSges(5,3);
t5 = -t21 * mrSges(6,1) + t19 * mrSges(6,2);
t85 = mrSges(5,1) * t75 + mrSges(5,2) * t76;
t84 = -Ifges(5,5) * t76 + Ifges(5,6) * t75;
t59 = t76 * t65;
t27 = -t76 * t105 + t88 * t78 + t59;
t30 = -t75 * t105 + t41;
t6 = t27 * t79 - t30 * t77;
t7 = t27 * t77 + t30 * t79;
t66 = t98 * t75;
t68 = t98 * t76;
t34 = -t66 * t79 - t68 * t77;
t35 = -t66 * t77 + t68 * t79;
t45 = t62 * t80;
t47 = t83 * t80;
t71 = -pkin(4) * t76 - pkin(3);
t64 = mrSges(5,1) * t78 - t76 * t99;
t63 = -mrSges(5,2) * t78 - t75 * t99;
t60 = t87 * t80;
t57 = (mrSges(5,1) * t80 + mrSges(5,3) * t101) * qJD(3);
t56 = (mrSges(5,3) * t75 * t78 - mrSges(5,2) * t80) * qJD(3);
t55 = t85 * t80;
t52 = t87 * t93;
t49 = t85 * t93;
t46 = t83 * t78;
t44 = t62 * t78;
t43 = t76 * t48;
t40 = -t75 * t100 + t59;
t39 = (t80 * Ifges(5,5) + (-t76 * Ifges(5,1) + t104) * t78) * qJD(3);
t38 = (t80 * Ifges(5,6) + (t102 - t103) * t78) * qJD(3);
t37 = mrSges(6,1) * t78 + mrSges(6,3) * t47;
t36 = -mrSges(6,2) * t78 - mrSges(6,3) * t45;
t33 = Ifges(6,1) * t62 - Ifges(6,4) * t83;
t32 = Ifges(6,4) * t62 - Ifges(6,2) * t83;
t31 = mrSges(6,1) * t83 + mrSges(6,2) * t62;
t28 = -t75 * t90 + t43;
t26 = -Ifges(6,1) * t112 - Ifges(6,4) * t54;
t25 = -Ifges(6,4) * t112 - Ifges(6,2) * t54;
t24 = mrSges(6,1) * t54 - mrSges(6,2) * t112;
t23 = mrSges(6,1) * t45 - mrSges(6,2) * t47;
t22 = pkin(7) * t75 * t93 + t29;
t20 = -qJD(3) * t45 + t112 * t78;
t18 = -qJD(3) * t47 - t78 * t54;
t14 = -Ifges(6,1) * t47 - Ifges(6,4) * t45 + Ifges(6,5) * t78;
t13 = -Ifges(6,4) * t47 - Ifges(6,2) * t45 + Ifges(6,6) * t78;
t12 = -t62 * qJD(4) - t35 * qJD(5);
t11 = -t83 * qJD(4) + t34 * qJD(5);
t10 = t43 + (pkin(7) * t101 + t80 * t88) * qJD(3);
t9 = -mrSges(6,2) * t92 + mrSges(6,3) * t21;
t8 = mrSges(6,1) * t92 - mrSges(6,3) * t19;
t4 = Ifges(6,1) * t19 + Ifges(6,4) * t21 + Ifges(6,5) * t92;
t3 = Ifges(6,4) * t19 + Ifges(6,2) * t21 + Ifges(6,6) * t92;
t2 = -qJD(5) * t7 + t10 * t79 - t22 * t77;
t1 = qJD(5) * t6 + t10 * t77 + t22 * t79;
t15 = [t78 * t91 + 0.2e1 * t29 * t63 + 0.2e1 * t28 * t64 - t47 * t4 - 0.2e1 * t52 * t23 + 0.2e1 * t41 * t56 + 0.2e1 * t40 * t57 + 0.2e1 * t60 * t5 + 0.2e1 * t1 * t36 + 0.2e1 * t2 * t37 - t45 * t3 + t19 * t14 + t21 * t13 + 0.2e1 * t6 * t8 + 0.2e1 * t7 * t9 + 0.2e1 * m(5) * (t40 * t28 + t41 * t29) + (t1 * t7 + t2 * t6 - t52 * t60) * t111 + (t114 * t49 - t75 * t38 + t76 * t39) * t80 + (t78 * t110 + 0.2e1 * t80 * mrSges(4,2) + (2 * mrSges(3,3)) + 0.2e1 * (m(3) + m(4)) * qJ(2)) * qJD(2) + ((qJ(2) * t110 - Ifges(6,5) * t47 - Ifges(6,6) * t45 + (-(2 * Ifges(4,4)) - t84) * t80) * t80 + (-0.2e1 * qJ(2) * mrSges(4,2) + t55 * t114 + 0.2e1 * (Ifges(4,4) + t84) * t78 + ((2 * Ifges(5,3)) - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + Ifges(6,3) - 0.2e1 * m(5) * (t81 ^ 2) - Ifges(5,1) * t74 + (-t102 + 0.2e1 * t103) * t75) * t80) * t78) * qJD(3); t18 * t36 + t20 * t37 - t44 * t8 - t46 * t9 + (t49 - t5) * t80 + (t76 * t56 - t75 * t57) * t78 + ((t76 * t63 - t75 * t64) * t80 + (t23 + t55) * t78) * qJD(3) + m(6) * (-t1 * t46 + t18 * t7 - t2 * t44 + t20 * t6 + t52 * t80 + t60 * t93) + m(5) * ((-t28 * t75 + t29 * t76) * t78 + (-t40 * t75 + t41 * t76 - 0.2e1 * t100) * t92); (-t18 * t46 - t20 * t44) * t111 + 0.4e1 * (m(5) * (-0.1e1 + t95) / 0.2e1 - m(6) / 0.2e1) * t78 * t92; t78 * t97 / 0.2e1 + t71 * t5 + t3 * t109 + t4 * t108 - t47 * t26 / 0.2e1 + pkin(3) * t49 - t52 * t31 - t112 * t14 / 0.2e1 - t54 * t13 / 0.2e1 + t60 * t24 + t35 * t9 + t11 * t36 + t12 * t37 - t45 * t25 / 0.2e1 + t21 * t32 / 0.2e1 + t19 * t33 / 0.2e1 + t34 * t8 + m(6) * (t1 * t35 + t11 * t7 + t12 * t6 + t2 * t34 - t52 * t71) + (t38 / 0.2e1 + qJD(4) * t63 + qJ(4) * t56 + t29 * mrSges(5,3) + m(5) * (qJ(4) * t29 + qJD(4) * t41)) * t76 + (t39 / 0.2e1 - qJD(4) * t64 - qJ(4) * t57 - t28 * mrSges(5,3) + m(5) * (-qJ(4) * t28 - qJD(4) * t40)) * t75 + (-t1 * t83 + t112 * t6 - t2 * t62 - t7 * t54) * mrSges(6,3) + ((Ifges(5,5) * t107 + Ifges(5,6) * t76 / 0.2e1 + Ifges(6,5) * t108 + Ifges(6,6) * t109 - Ifges(4,6) - t81 * mrSges(4,2)) * t80 + (-Ifges(4,5) - t76 * (Ifges(5,1) * t75 + t103) / 0.2e1 + (Ifges(5,2) * t76 + t104) * t107 + (-m(5) * pkin(3) + t96) * t81) * t78) * qJD(3); -t80 * t24 + m(6) * (-t11 * t46 - t12 * t44 + t18 * t35 + t20 * t34) + t78 * t113 + ((-mrSges(4,2) + t86) * t80 + m(5) * (t94 * t95 - t106) + (m(6) * t71 + t31 + t96) * t78) * qJD(3) + (-t112 * t44 - t18 * t83 - t20 * t62 + t46 * t54) * mrSges(6,3); -t112 * t33 + t62 * t26 - t54 * t32 - t83 * t25 + 0.2e1 * t71 * t24 + (t11 * t35 + t12 * t34) * t111 + 0.2e1 * qJ(4) * t113 + 0.2e1 * t86 * qJD(4) + 0.2e1 * (-t11 * t83 + t112 * t34 - t12 * t62 - t35 * t54) * mrSges(6,3); -m(6) * t52 + (m(5) * t81 - t85) * t93 + t5; (m(5) + m(6)) * t93; t24; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t91; mrSges(6,1) * t20 - mrSges(6,2) * t18; mrSges(6,1) * t12 - mrSges(6,2) * t11 + t97; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;

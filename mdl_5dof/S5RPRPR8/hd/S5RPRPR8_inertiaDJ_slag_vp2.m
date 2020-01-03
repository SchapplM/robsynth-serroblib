% Calculate time derivative of joint inertia matrix for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:19
% EndTime: 2019-12-31 18:21:22
% DurationCPUTime: 1.17s
% Computational Cost: add. (1136->223), mult. (2741->358), div. (0->0), fcn. (2247->8), ass. (0->99)
t73 = sin(pkin(9));
t74 = cos(pkin(9));
t116 = mrSges(5,1) * t73 + mrSges(5,2) * t74;
t85 = -cos(pkin(8)) * pkin(1) - pkin(2);
t115 = 0.2e1 * t85;
t72 = t74 ^ 2;
t92 = t73 ^ 2 + t72;
t76 = sin(qJ(5));
t78 = cos(qJ(5));
t80 = t73 * t76 - t74 * t78;
t51 = t80 * qJD(5);
t113 = 2 * m(5);
t112 = 2 * m(6);
t68 = sin(pkin(8)) * pkin(1) + pkin(6);
t111 = 0.2e1 * t68;
t110 = -t80 / 0.2e1;
t59 = t73 * t78 + t74 * t76;
t109 = t59 / 0.2e1;
t108 = t74 / 0.2e1;
t77 = sin(qJ(3));
t107 = pkin(3) * t77;
t79 = cos(qJ(3));
t89 = qJD(3) * t79;
t46 = t116 * t89;
t52 = t59 * qJD(5);
t19 = -t77 * t52 - t80 * t89;
t20 = t51 * t77 - t59 * t89;
t7 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t106 = -t46 - t7;
t103 = Ifges(5,4) * t73;
t102 = Ifges(5,4) * t74;
t101 = t73 * Ifges(5,2);
t100 = t73 * t77;
t99 = t73 * t79;
t98 = t74 * t77;
t97 = t74 * t79;
t96 = t79 * mrSges(5,3);
t95 = pkin(7) + qJ(4);
t94 = -Ifges(6,5) * t51 - Ifges(6,6) * t52;
t88 = t77 * qJD(4);
t91 = qJ(4) * t79;
t50 = -t88 + (-t91 + t107) * qJD(3);
t90 = qJD(3) * t77;
t86 = t68 * t90;
t26 = t74 * t50 + t73 * t86;
t56 = -pkin(3) * t79 - t77 * qJ(4) + t85;
t29 = t73 * t56 + t68 * t97;
t93 = -mrSges(5,1) * t74 + mrSges(5,2) * t73 - mrSges(4,1);
t87 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t90;
t84 = m(5) * t92;
t83 = pkin(4) * t73 + t68;
t82 = t92 * mrSges(5,3);
t81 = -Ifges(5,5) * t74 + Ifges(5,6) * t73;
t45 = t74 * t56;
t18 = -pkin(7) * t98 + t45 + (-t68 * t73 - pkin(4)) * t79;
t22 = -pkin(7) * t100 + t29;
t5 = t18 * t78 - t22 * t76;
t6 = t18 * t76 + t22 * t78;
t63 = t95 * t73;
t65 = t95 * t74;
t33 = -t63 * t78 - t65 * t76;
t34 = -t63 * t76 + t65 * t78;
t69 = -pkin(4) * t74 - pkin(3);
t61 = -mrSges(5,1) * t79 - mrSges(5,3) * t98;
t60 = mrSges(5,2) * t79 - mrSges(5,3) * t100;
t55 = (mrSges(5,1) * t77 - t74 * t96) * qJD(3);
t54 = (-mrSges(5,2) * t77 - t73 * t96) * qJD(3);
t53 = t116 * t77;
t49 = t83 * t77;
t43 = t80 * t77;
t42 = t59 * t77;
t41 = t83 * t89;
t39 = t73 * t50;
t38 = (t77 * Ifges(5,5) + (t74 * Ifges(5,1) - t103) * t79) * qJD(3);
t37 = (t77 * Ifges(5,6) + (-t101 + t102) * t79) * qJD(3);
t36 = -mrSges(6,1) * t79 + t43 * mrSges(6,3);
t35 = mrSges(6,2) * t79 - t42 * mrSges(6,3);
t32 = Ifges(6,1) * t59 - Ifges(6,4) * t80;
t31 = Ifges(6,4) * t59 - Ifges(6,2) * t80;
t30 = mrSges(6,1) * t80 + mrSges(6,2) * t59;
t28 = -t68 * t99 + t45;
t27 = -t74 * t86 + t39;
t25 = -Ifges(6,1) * t51 - Ifges(6,4) * t52;
t24 = -Ifges(6,4) * t51 - Ifges(6,2) * t52;
t23 = mrSges(6,1) * t52 - mrSges(6,2) * t51;
t21 = mrSges(6,1) * t42 - mrSges(6,2) * t43;
t15 = t39 + (-pkin(7) * t99 - t68 * t98) * qJD(3);
t14 = -Ifges(6,1) * t43 - Ifges(6,4) * t42 - Ifges(6,5) * t79;
t13 = -Ifges(6,4) * t43 - Ifges(6,2) * t42 - Ifges(6,6) * t79;
t12 = -t59 * qJD(4) - t34 * qJD(5);
t11 = -t80 * qJD(4) + t33 * qJD(5);
t10 = (pkin(4) * t77 - pkin(7) * t97) * qJD(3) + t26;
t9 = -mrSges(6,2) * t90 + mrSges(6,3) * t20;
t8 = mrSges(6,1) * t90 - mrSges(6,3) * t19;
t4 = Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t90;
t3 = Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t90;
t2 = -t6 * qJD(5) + t10 * t78 - t15 * t76;
t1 = t5 * qJD(5) + t10 * t76 + t15 * t78;
t16 = [-t79 * t87 + 0.2e1 * t29 * t54 + 0.2e1 * t28 * t55 + 0.2e1 * t27 * t60 + 0.2e1 * t26 * t61 - t42 * t3 - t43 * t4 + 0.2e1 * t49 * t7 + 0.2e1 * t1 * t35 + 0.2e1 * t2 * t36 + 0.2e1 * t41 * t21 + t19 * t14 + t20 * t13 + 0.2e1 * t5 * t8 + 0.2e1 * t6 * t9 + (t28 * t26 + t29 * t27) * t113 + (t1 * t6 + t2 * t5 + t41 * t49) * t112 + (t46 * t111 - t73 * t37 + t74 * t38) * t77 + ((-Ifges(6,5) * t43 - Ifges(6,6) * t42 + mrSges(4,1) * t115 + (-(2 * Ifges(4,4)) - t81) * t77) * t77 + (t53 * t111 + mrSges(4,2) * t115 + 0.2e1 * (Ifges(4,4) + t81) * t79 + (Ifges(5,1) * t72 - (2 * Ifges(5,3)) + t68 ^ 2 * t113 - Ifges(6,3) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + (t101 - 0.2e1 * t102) * t73) * t77) * t79) * qJD(3); t19 * t35 + t20 * t36 - t42 * t8 - t43 * t9 + t106 * t79 + (t74 * t54 - t73 * t55) * t77 + ((t60 * t74 - t61 * t73) * t79 + (t21 + t53) * t77) * qJD(3) + m(6) * (-t1 * t43 + t6 * t19 - t2 * t42 + t5 * t20 - t41 * t79 + t49 * t90) + m(5) * ((-t26 * t73 + t27 * t74) * t77 + (t68 * t77 ^ 2 + (-t28 * t73 + t29 * t74 - t68 * t79) * t79) * qJD(3)); (-t43 * t19 - t42 * t20) * t112 + 0.4e1 * (-m(6) / 0.2e1 + m(5) * (-0.1e1 + t92) / 0.2e1) * t77 * t89; m(6) * (t1 * t34 + t11 * t6 + t12 * t5 + t2 * t33 + t41 * t69) - t79 * t94 / 0.2e1 + t69 * t7 + t3 * t110 + t4 * t109 - t42 * t24 / 0.2e1 - t43 * t25 / 0.2e1 - pkin(3) * t46 + t49 * t23 - t51 * t14 / 0.2e1 - t52 * t13 / 0.2e1 + t20 * t31 / 0.2e1 + t19 * t32 / 0.2e1 + t33 * t8 + t34 * t9 + t11 * t35 + t12 * t36 + t41 * t30 + (m(5) * (qJ(4) * t27 + qJD(4) * t29) + qJD(4) * t60 + qJ(4) * t54 + t27 * mrSges(5,3) + t37 / 0.2e1) * t74 + (m(5) * (-qJ(4) * t26 - qJD(4) * t28) - qJD(4) * t61 - qJ(4) * t55 - t26 * mrSges(5,3) + t38 / 0.2e1) * t73 + (-t1 * t80 - t2 * t59 + t5 * t51 - t6 * t52) * mrSges(6,3) + ((t68 * mrSges(4,2) + Ifges(5,5) * t73 / 0.2e1 + Ifges(5,6) * t108 + Ifges(6,5) * t109 + Ifges(6,6) * t110 - Ifges(4,6)) * t77 + ((Ifges(5,1) * t73 + t102) * t108 - t73 * (Ifges(5,2) * t74 + t103) / 0.2e1 + Ifges(4,5) + (-m(5) * pkin(3) + t93) * t68) * t79) * qJD(3); -t79 * t23 + m(6) * (-t11 * t43 - t12 * t42 + t19 * t34 + t20 * t33) + t84 * t88 + ((-mrSges(4,2) + t82) * t79 + m(5) * (t91 * t92 - t107) + (m(6) * t69 + t30 + t93) * t77) * qJD(3) + (-t19 * t80 - t20 * t59 - t42 * t51 + t43 * t52) * mrSges(6,3); (t11 * t34 + t12 * t33) * t112 - t51 * t32 + t59 * t25 - t52 * t31 - t80 * t24 + 0.2e1 * t69 * t23 + 0.2e1 * (-t11 * t80 - t12 * t59 + t33 * t51 - t34 * t52) * mrSges(6,3) + 0.2e1 * (qJ(4) * t84 + t82) * qJD(4); m(5) * t68 * t89 + m(6) * t41 - t106; (m(5) + m(6)) * t90; t23; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t87; -t7; mrSges(6,1) * t12 - mrSges(6,2) * t11 + t94; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR13_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:58
% EndTime: 2019-12-31 18:32:00
% DurationCPUTime: 0.86s
% Computational Cost: add. (1137->170), mult. (2537->250), div. (0->0), fcn. (2295->6), ass. (0->84)
t71 = 2 * mrSges(5,1) + 2 * mrSges(4,3);
t54 = cos(pkin(8));
t94 = cos(qJ(3));
t75 = t94 * t54;
t53 = sin(pkin(8));
t56 = sin(qJ(3));
t89 = t53 * t56;
t37 = -t75 + t89;
t57 = cos(qJ(5));
t80 = qJD(5) * t57;
t38 = t94 * t53 + t56 * t54;
t34 = t38 * qJD(3);
t55 = sin(qJ(5));
t88 = t55 * t34;
t61 = t37 * t80 + t88;
t83 = pkin(6) + qJ(2);
t42 = t83 * t53;
t43 = t83 * t54;
t26 = -t56 * t42 + t94 * t43;
t18 = t38 * qJD(2) + t26 * qJD(3);
t74 = qJD(3) * t94;
t33 = qJD(3) * t89 - t54 * t74;
t10 = -t33 * pkin(4) + t18;
t76 = -pkin(2) * t54 - pkin(1);
t62 = -t38 * qJ(4) + t76;
t97 = pkin(3) + pkin(7);
t14 = t97 * t37 + t62;
t25 = t94 * t42 + t43 * t56;
t19 = t38 * pkin(4) + t25;
t3 = -t14 * t55 + t19 * t57;
t64 = qJ(4) * t33 - qJD(4) * t38;
t8 = t97 * t34 + t64;
t1 = qJD(5) * t3 + t10 * t55 + t57 * t8;
t4 = t14 * t57 + t19 * t55;
t2 = -qJD(5) * t4 + t10 * t57 - t55 * t8;
t103 = t1 * t55 + t2 * t57;
t17 = (qJD(2) * t53 + qJD(3) * t43) * t56 - qJD(2) * t75 + t42 * t74;
t13 = pkin(3) * t34 + t64;
t102 = -0.2e1 * t13;
t101 = -t55 / 0.2e1;
t100 = t55 / 0.2e1;
t99 = -t57 / 0.2e1;
t98 = t57 / 0.2e1;
t93 = mrSges(6,3) * t37;
t92 = Ifges(6,4) * t55;
t91 = Ifges(6,4) * t57;
t90 = Ifges(6,6) * t38;
t46 = Ifges(6,1) * t57 - t92;
t87 = t55 * t46;
t86 = t57 * t34;
t45 = -Ifges(6,2) * t55 + t91;
t85 = t57 * t45;
t84 = mrSges(5,2) - mrSges(4,1);
t81 = qJD(5) * t55;
t79 = 0.2e1 * t34;
t78 = t37 * t81;
t72 = t61 * Ifges(6,5) + Ifges(6,6) * t86 - Ifges(6,3) * t33;
t70 = t3 * t55 - t4 * t57;
t44 = mrSges(6,1) * t55 + mrSges(6,2) * t57;
t69 = Ifges(6,1) * t55 + t91;
t68 = Ifges(6,2) * t57 + t92;
t67 = -Ifges(6,5) * t55 - Ifges(6,6) * t57;
t66 = -t17 * t26 + t18 * t25;
t23 = t38 * mrSges(6,1) - t55 * t93;
t24 = -t38 * mrSges(6,2) + t57 * t93;
t65 = -t55 * t23 + t57 * t24;
t60 = t78 - t86;
t41 = t69 * qJD(5);
t40 = t68 * qJD(5);
t39 = -mrSges(6,1) * t80 + mrSges(6,2) * t81;
t31 = t33 * mrSges(4,2);
t30 = t33 * mrSges(5,3);
t22 = (-mrSges(6,1) * t57 + mrSges(6,2) * t55) * t37;
t21 = t37 * pkin(3) + t62;
t20 = -t37 * pkin(4) + t26;
t16 = Ifges(6,5) * t38 + t69 * t37;
t15 = t68 * t37 + t90;
t12 = -t33 * mrSges(6,1) - t61 * mrSges(6,3);
t11 = t33 * mrSges(6,2) - t60 * mrSges(6,3);
t9 = -t34 * pkin(4) - t17;
t7 = t60 * mrSges(6,1) + t61 * mrSges(6,2);
t6 = t61 * Ifges(6,1) - t60 * Ifges(6,4) - Ifges(6,5) * t33;
t5 = t61 * Ifges(6,4) - t60 * Ifges(6,2) - Ifges(6,6) * t33;
t27 = [t15 * t86 + t16 * t88 + 0.2e1 * t76 * (t34 * mrSges(4,1) - t31) + 0.2e1 * t21 * (-t34 * mrSges(5,2) + t30) + 0.2e1 * t20 * t7 + 0.2e1 * t9 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t1 * t24 + 0.2e1 * t4 * t11 + 0.2e1 * t3 * t12 + 0.2e1 * m(6) * (t1 * t4 + t2 * t3 + t20 * t9) + 0.2e1 * m(4) * t66 + 0.2e1 * m(5) * (t13 * t21 + t66) + (mrSges(5,3) * t102 + (-Ifges(4,4) - Ifges(5,6)) * t79 + t18 * t71 + (-(2 * Ifges(4,1)) - (2 * Ifges(5,2)) - Ifges(6,3)) * t33 + t72) * t38 + (mrSges(5,2) * t102 + t57 * t5 + t55 * t6 + (Ifges(5,3) + Ifges(4,2)) * t79 + t17 * t71 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t67) * t33 + (t57 * t16 + (-t15 - t90) * t55) * qJD(5)) * t37 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t53 ^ 2 + t54 ^ 2) * qJD(2) + (-t25 * t33 - t26 * t34) * t71; t57 * t11 - t55 * t12 + t30 - t31 - t84 * t34 + (-t57 * t23 - t55 * t24) * qJD(5) + m(6) * (t1 * t57 - t2 * t55 + (-t3 * t57 - t4 * t55) * qJD(5)) + m(5) * t13; 0; qJ(4) * t7 + qJD(4) * t22 - t20 * t39 + t9 * t44 + t84 * t18 + (-mrSges(5,3) + mrSges(4,2)) * t17 + (-t2 * mrSges(6,3) - t97 * t12 + t6 / 0.2e1) * t57 + (-t97 * t11 - t1 * mrSges(6,3) - t5 / 0.2e1) * t55 + m(6) * (qJ(4) * t9 + qJD(4) * t20 - t103 * t97) + m(5) * (-pkin(3) * t18 - qJ(4) * t17 + qJD(4) * t26) + (-qJD(4) * mrSges(5,1) - t41 * t100 - t40 * t98) * t37 + (pkin(3) * mrSges(5,1) + Ifges(6,5) * t99 + Ifges(6,6) * t100 + Ifges(5,4) - Ifges(4,5)) * t33 + (t85 / 0.2e1 + t87 / 0.2e1 - qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6)) * t34 + (t15 * t99 + t16 * t101 + t38 * t67 / 0.2e1 + (t45 * t101 + t46 * t98) * t37 + t70 * mrSges(6,3) - (-m(6) * t70 + t65) * t97) * qJD(5); 0; -0.2e1 * qJ(4) * t39 + t40 * t55 - t41 * t57 + (-t85 - t87) * qJD(5) + 0.2e1 * (mrSges(5,3) + t44 + (m(5) + m(6)) * qJ(4)) * qJD(4); -t33 * mrSges(5,1) + t55 * t11 + t57 * t12 + t65 * qJD(5) + m(6) * (-t70 * qJD(5) + t103) + m(5) * t18; 0; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) - Ifges(6,6) * t78 + t72; t39; ((mrSges(6,2) * t97 - Ifges(6,6)) * t57 + (mrSges(6,1) * t97 - Ifges(6,5)) * t55) * qJD(5); -t44 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t27(1), t27(2), t27(4), t27(7), t27(11); t27(2), t27(3), t27(5), t27(8), t27(12); t27(4), t27(5), t27(6), t27(9), t27(13); t27(7), t27(8), t27(9), t27(10), t27(14); t27(11), t27(12), t27(13), t27(14), t27(15);];
Mq = res;

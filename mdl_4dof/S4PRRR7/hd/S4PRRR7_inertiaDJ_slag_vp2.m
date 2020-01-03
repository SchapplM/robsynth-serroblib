% Calculate time derivative of joint inertia matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:04
% EndTime: 2019-12-31 16:36:06
% DurationCPUTime: 0.64s
% Computational Cost: add. (424->145), mult. (1248->248), div. (0->0), fcn. (981->8), ass. (0->78)
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t29 = -t41 * mrSges(5,1) + t38 * mrSges(5,2);
t82 = -m(5) * pkin(3) - mrSges(4,1) + t29;
t81 = 0.2e1 * m(5);
t80 = 2 * pkin(6);
t72 = Ifges(5,4) * t38;
t30 = Ifges(5,2) * t41 + t72;
t79 = -t30 / 0.2e1;
t78 = -t38 / 0.2e1;
t42 = cos(qJ(3));
t77 = pkin(6) * t42;
t37 = cos(pkin(4));
t39 = sin(qJ(3));
t36 = sin(pkin(4));
t40 = sin(qJ(2));
t68 = t36 * t40;
t18 = -t37 * t42 + t39 * t68;
t19 = t37 * t39 + t42 * t68;
t43 = cos(qJ(2));
t59 = t36 * qJD(2);
t58 = t43 * t59;
t8 = t19 * qJD(3) + t39 * t58;
t76 = t18 * t8;
t75 = t8 * t39;
t9 = -t18 * qJD(3) + t42 * t58;
t74 = t9 * t42;
t73 = mrSges(5,3) * t39;
t71 = Ifges(5,4) * t41;
t70 = Ifges(5,6) * t38;
t69 = Ifges(5,6) * t42;
t67 = t36 * t43;
t28 = -t42 * pkin(3) - t39 * pkin(7) - pkin(2);
t66 = t38 * t28;
t63 = qJD(3) * t42;
t55 = t41 * t63;
t64 = qJD(3) * t39;
t65 = -Ifges(5,5) * t55 - Ifges(5,3) * t64;
t62 = qJD(4) * t39;
t61 = qJD(4) * t41;
t60 = qJD(4) * t42;
t57 = t38 * t62;
t56 = t40 * t59;
t54 = (2 * Ifges(4,4)) + t70;
t14 = t41 * t28 - t38 * t77;
t27 = (pkin(3) * t39 - pkin(7) * t42) * qJD(3);
t3 = t28 * t61 + t38 * t27 + (-t38 * t60 - t41 * t64) * pkin(6);
t53 = -t14 * qJD(4) + t3;
t52 = mrSges(5,1) * t38 + mrSges(5,2) * t41;
t51 = Ifges(5,1) * t41 - t72;
t31 = Ifges(5,1) * t38 + t71;
t50 = -Ifges(5,2) * t38 + t71;
t49 = Ifges(5,5) * t38 + Ifges(5,6) * t41;
t10 = -t19 * t38 - t41 * t67;
t48 = -t19 * t41 + t38 * t67;
t47 = t18 * t63 + t75;
t46 = t55 - t57;
t45 = t38 * t63 + t39 * t61;
t35 = Ifges(5,5) * t61;
t26 = -t42 * mrSges(5,1) - t41 * t73;
t25 = t42 * mrSges(5,2) - t38 * t73;
t24 = t51 * qJD(4);
t23 = t50 * qJD(4);
t22 = (mrSges(4,1) * t39 + mrSges(4,2) * t42) * qJD(3);
t21 = t52 * qJD(4);
t20 = t52 * t39;
t17 = -Ifges(5,5) * t42 + t51 * t39;
t16 = t50 * t39 - t69;
t15 = t41 * t77 + t66;
t13 = -mrSges(5,2) * t64 - t45 * mrSges(5,3);
t12 = mrSges(5,1) * t64 - t46 * mrSges(5,3);
t7 = t45 * mrSges(5,1) + t46 * mrSges(5,2);
t6 = -t31 * t62 + (Ifges(5,5) * t39 + t51 * t42) * qJD(3);
t5 = -t30 * t62 + (Ifges(5,6) * t39 + t50 * t42) * qJD(3);
t4 = -qJD(4) * t66 + t41 * t27 + (t38 * t64 - t41 * t60) * pkin(6);
t2 = t10 * qJD(4) + t38 * t56 + t9 * t41;
t1 = t48 * qJD(4) - t9 * t38 + t41 * t56;
t11 = [0.2e1 * m(4) * (-t36 ^ 2 * t40 * qJD(2) * t43 + t19 * t9 + t76) + 0.2e1 * m(5) * (t10 * t1 - t2 * t48 + t76); t1 * t26 + t10 * t12 - t48 * t13 + t18 * t7 + t2 * t25 + t8 * t20 + (-t43 * t22 + (-t43 * mrSges(3,2) + (-t42 * mrSges(4,1) + t39 * mrSges(4,2) - mrSges(3,1)) * t40) * qJD(2)) * t36 - m(4) * pkin(2) * t56 + m(5) * (t14 * t1 + t4 * t10 + t15 * t2 - t3 * t48) + (m(4) * (-t19 * t64 + t47 + t74) / 0.2e1 + m(5) * t47 / 0.2e1) * t80 + (t75 + t74 + (t18 * t42 - t19 * t39) * qJD(3)) * mrSges(4,3); (t14 * t4 + t15 * t3) * t81 + 0.2e1 * t3 * t25 + 0.2e1 * t15 * t13 + 0.2e1 * t4 * t26 + 0.2e1 * t14 * t12 - 0.2e1 * pkin(2) * t22 + ((-t38 * t16 + t41 * t17 + t20 * t80 + t54 * t42) * qJD(3) + t65) * t42 + (t7 * t80 - t38 * t5 + t41 * t6 + (-t41 * t16 - t38 * t17 + t42 * t49) * qJD(4) + ((Ifges(5,5) * t41 - t54) * t39 + ((pkin(6) ^ 2) * t81 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3)) * t42) * qJD(3)) * t39; -t9 * mrSges(4,2) + t18 * t21 + (m(5) * pkin(7) + mrSges(5,3)) * (-t1 * t38 + t2 * t41 + (-t10 * t41 + t38 * t48) * qJD(4)) + t82 * t8; -pkin(3) * t7 + (-t35 / 0.2e1 + (t82 * pkin(6) + Ifges(4,5)) * qJD(3)) * t42 + (t6 / 0.2e1 + t63 * t79 - t4 * mrSges(5,3) + (t69 / 0.2e1 - t15 * mrSges(5,3) - t16 / 0.2e1) * qJD(4) + (-qJD(4) * t25 + m(5) * (-t15 * qJD(4) - t4) - t12) * pkin(7)) * t38 + (t5 / 0.2e1 + t31 * t63 / 0.2e1 + qJD(4) * t17 / 0.2e1 + t53 * mrSges(5,3) + (m(5) * t53 - qJD(4) * t26 + t13) * pkin(7)) * t41 + (t23 * t78 + t41 * t24 / 0.2e1 + (t31 * t78 + t41 * t79) * qJD(4) + pkin(6) * t21 + (-Ifges(4,6) + t49 / 0.2e1 + pkin(6) * mrSges(4,2)) * qJD(3)) * t39; -0.2e1 * pkin(3) * t21 + t41 * t23 + t38 * t24 + (-t30 * t38 + t31 * t41) * qJD(4); t1 * mrSges(5,1) - t2 * mrSges(5,2); t4 * mrSges(5,1) - t3 * mrSges(5,2) - Ifges(5,5) * t57 - t45 * Ifges(5,6) - t65; t35 + (t29 * pkin(7) - t70) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1), t11(2), t11(4), t11(7); t11(2), t11(3), t11(5), t11(8); t11(4), t11(5), t11(6), t11(9); t11(7), t11(8), t11(9), t11(10);];
Mq = res;

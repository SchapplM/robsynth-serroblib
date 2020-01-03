% Calculate time derivative of joint inertia matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:58
% EndTime: 2019-12-31 17:06:00
% DurationCPUTime: 0.67s
% Computational Cost: add. (741->145), mult. (1722->226), div. (0->0), fcn. (1454->6), ass. (0->73)
t81 = -2 * mrSges(4,3);
t46 = cos(qJ(2));
t66 = -qJ(3) - pkin(5);
t33 = t66 * t46;
t41 = sin(pkin(7));
t42 = cos(pkin(7));
t44 = sin(qJ(2));
t57 = t66 * t44;
t18 = -t33 * t41 - t42 * t57;
t79 = 0.2e1 * t18;
t28 = t41 * t46 + t42 * t44;
t78 = -t28 / 0.2e1;
t45 = cos(qJ(4));
t43 = sin(qJ(4));
t75 = Ifges(5,4) * t43;
t34 = Ifges(5,2) * t45 + t75;
t77 = -t34 / 0.2e1;
t76 = mrSges(5,3) * t28;
t74 = Ifges(5,4) * t45;
t25 = t28 * qJD(2);
t73 = Ifges(5,5) * t25;
t72 = Ifges(5,6) * t25;
t27 = t41 * t44 - t42 * t46;
t71 = Ifges(5,6) * t27;
t70 = Ifges(5,6) * t43;
t53 = qJD(2) * t66;
t24 = qJD(3) * t46 + t44 * t53;
t47 = -t44 * qJD(3) + t46 * t53;
t12 = t24 * t41 - t42 * t47;
t69 = t12 * t18;
t26 = t27 * qJD(2);
t67 = t45 * t26;
t65 = -Ifges(5,5) * t67 + Ifges(5,3) * t25;
t39 = -pkin(2) * t46 - pkin(1);
t15 = t27 * pkin(3) - t28 * pkin(6) + t39;
t19 = -t42 * t33 + t41 * t57;
t6 = t15 * t45 - t19 * t43;
t64 = qJD(4) * t6;
t7 = t15 * t43 + t19 * t45;
t63 = qJD(4) * t7;
t62 = qJD(4) * t43;
t61 = qJD(4) * t45;
t59 = pkin(2) * qJD(2) * t44;
t58 = t28 * t62;
t56 = t25 * mrSges(4,1) - t26 * mrSges(4,2);
t55 = -(2 * Ifges(4,4)) - t70;
t54 = 0.2e1 * t59;
t52 = mrSges(5,1) * t43 + mrSges(5,2) * t45;
t51 = Ifges(5,1) * t45 - t75;
t50 = -Ifges(5,2) * t43 + t74;
t49 = -t26 * t43 + t28 * t61;
t48 = t58 + t67;
t40 = Ifges(5,5) * t61;
t38 = -pkin(2) * t42 - pkin(3);
t37 = pkin(2) * t41 + pkin(6);
t35 = Ifges(5,1) * t43 + t74;
t32 = t51 * qJD(4);
t31 = t50 * qJD(4);
t30 = t52 * qJD(4);
t17 = mrSges(5,1) * t27 - t45 * t76;
t16 = -mrSges(5,2) * t27 - t43 * t76;
t14 = pkin(3) * t25 + pkin(6) * t26 + t59;
t13 = t42 * t24 + t41 * t47;
t11 = Ifges(5,5) * t27 + t51 * t28;
t10 = t50 * t28 + t71;
t9 = -mrSges(5,2) * t25 - t49 * mrSges(5,3);
t8 = mrSges(5,1) * t25 + t48 * mrSges(5,3);
t5 = t49 * mrSges(5,1) - t48 * mrSges(5,2);
t4 = -t48 * Ifges(5,1) - t49 * Ifges(5,4) + t73;
t3 = -t48 * Ifges(5,4) - t49 * Ifges(5,2) + t72;
t2 = -t13 * t43 + t14 * t45 - t63;
t1 = t13 * t45 + t14 * t43 + t64;
t20 = [t19 * t25 * t81 + 0.2e1 * t39 * t56 + 0.2e1 * t1 * t16 + 0.2e1 * t2 * t17 + t5 * t79 + 0.2e1 * t6 * t8 + 0.2e1 * t7 * t9 + 0.2e1 * m(4) * (t13 * t19 + t39 * t59 + t69) + 0.2e1 * m(5) * (t1 * t7 + t2 * t6 + t69) - (mrSges(4,3) * t79 - t43 * t10 + t45 * t11) * t26 + (mrSges(4,1) * t54 + t13 * t81 - t55 * t26 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t25 + t65) * t27 + (-t43 * t3 + t45 * t4 - 0.2e1 * Ifges(4,1) * t26 + mrSges(4,2) * t54 + (Ifges(5,5) * t45 + t55) * t25 + (t27 * (-Ifges(5,5) * t43 - Ifges(5,6) * t45) - t43 * t11 - t45 * t10) * qJD(4) + 0.2e1 * (mrSges(4,3) + t52) * t12) * t28 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t44 + mrSges(3,2) * t46) + (Ifges(3,1) - Ifges(3,2)) * t44 * t46 + (-t44 ^ 2 + t46 ^ 2) * Ifges(3,4)) * qJD(2); t38 * t5 + t27 * t40 / 0.2e1 + t18 * t30 - Ifges(4,6) * t25 - Ifges(4,5) * t26 - t13 * mrSges(4,2) + (m(5) * t38 - mrSges(4,1)) * t12 + (Ifges(3,5) * t46 - Ifges(3,6) * t44 + (-mrSges(3,1) * t46 + mrSges(3,2) * t44) * pkin(5)) * qJD(2) + (m(4) * (-t12 * t42 + t13 * t41) + (-t25 * t41 + t26 * t42) * mrSges(4,3)) * pkin(2) + (t31 * t78 - t26 * t77 - t2 * mrSges(5,3) + t73 / 0.2e1 + t12 * mrSges(5,2) + t4 / 0.2e1 + (-t10 / 0.2e1 - t71 / 0.2e1 + t35 * t78 - t7 * mrSges(5,3)) * qJD(4) + (-t8 + m(5) * (-t2 - t63) - qJD(4) * t16) * t37) * t43 + (t28 * t32 / 0.2e1 - t26 * t35 / 0.2e1 + t1 * mrSges(5,3) + t72 / 0.2e1 - t12 * mrSges(5,1) + t3 / 0.2e1 + (t11 / 0.2e1 + t28 * t77 - t6 * mrSges(5,3)) * qJD(4) + (m(5) * (t1 - t64) + t9 - qJD(4) * t17) * t37) * t45; 0.2e1 * t30 * t38 + t31 * t45 + t32 * t43 + (-t43 * t34 + t45 * t35) * qJD(4); m(5) * (t1 * t43 + t2 * t45 + (-t43 * t6 + t45 * t7) * qJD(4)) + t16 * t61 + t43 * t9 - t17 * t62 + t45 * t8 + m(4) * t59 + t56; 0; 0; mrSges(5,1) * t2 - mrSges(5,2) * t1 - Ifges(5,5) * t58 - t49 * Ifges(5,6) + t65; t40 + (-t70 + (-t45 * mrSges(5,1) + t43 * mrSges(5,2)) * t37) * qJD(4); -t30; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t20(1), t20(2), t20(4), t20(7); t20(2), t20(3), t20(5), t20(8); t20(4), t20(5), t20(6), t20(9); t20(7), t20(8), t20(9), t20(10);];
Mq = res;

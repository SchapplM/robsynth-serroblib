% Calculate time derivative of joint inertia matrix for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:15
% EndTime: 2019-12-31 16:50:16
% DurationCPUTime: 0.50s
% Computational Cost: add. (359->122), mult. (886->198), div. (0->0), fcn. (594->6), ass. (0->66)
t41 = -cos(pkin(7)) * pkin(1) - pkin(2);
t70 = 0.2e1 * t41;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t20 = -t32 * mrSges(5,1) + t30 * mrSges(5,2);
t69 = -mrSges(4,1) + t20;
t52 = t30 ^ 2 + t32 ^ 2;
t68 = 2 * m(5);
t24 = sin(pkin(7)) * pkin(1) + pkin(5);
t67 = 0.2e1 * t24;
t61 = Ifges(5,4) * t30;
t21 = Ifges(5,2) * t32 + t61;
t66 = -t21 / 0.2e1;
t65 = -t30 / 0.2e1;
t31 = sin(qJ(3));
t64 = pkin(3) * t31;
t33 = cos(qJ(3));
t63 = pkin(6) * t33;
t62 = mrSges(5,3) * t31;
t60 = Ifges(5,4) * t32;
t59 = Ifges(5,5) * t30;
t58 = Ifges(5,6) * t30;
t57 = Ifges(5,6) * t32;
t56 = Ifges(5,6) * t33;
t55 = t30 * t33;
t54 = t32 * t33;
t48 = qJD(3) * t33;
t42 = t32 * t48;
t49 = qJD(3) * t31;
t53 = -Ifges(5,5) * t42 - Ifges(5,3) * t49;
t13 = -t33 * pkin(3) - t31 * pkin(6) + t41;
t6 = t32 * t13 - t24 * t55;
t51 = qJD(4) * t6;
t7 = t30 * t13 + t24 * t54;
t50 = qJD(4) * t7;
t47 = qJD(4) * t30;
t46 = qJD(4) * t31;
t45 = qJD(4) * t32;
t44 = t24 * t49;
t43 = t30 * t46;
t40 = (2 * Ifges(4,4)) + t58;
t19 = (-t63 + t64) * qJD(3);
t1 = t30 * t19 - t32 * t44 + t51;
t39 = t1 - t51;
t38 = mrSges(5,1) * t30 + mrSges(5,2) * t32;
t37 = Ifges(5,1) * t32 - t61;
t22 = Ifges(5,1) * t30 + t60;
t36 = -Ifges(5,2) * t30 + t60;
t35 = t42 - t43;
t34 = t30 * t48 + t31 * t45;
t26 = Ifges(5,5) * t45;
t18 = -t33 * mrSges(5,1) - t32 * t62;
t17 = t33 * mrSges(5,2) - t30 * t62;
t16 = t37 * qJD(4);
t15 = t36 * qJD(4);
t14 = t38 * qJD(4);
t12 = t38 * t31;
t11 = -Ifges(5,5) * t33 + t37 * t31;
t10 = t36 * t31 - t56;
t9 = -mrSges(5,2) * t49 - t34 * mrSges(5,3);
t8 = mrSges(5,1) * t49 - t35 * mrSges(5,3);
t5 = t34 * mrSges(5,1) + t35 * mrSges(5,2);
t4 = -t22 * t46 + (Ifges(5,5) * t31 + t37 * t33) * qJD(3);
t3 = -t21 * t46 + (Ifges(5,6) * t31 + t36 * t33) * qJD(3);
t2 = t32 * t19 + t30 * t44 - t50;
t23 = [(t7 * t1 + t6 * t2) * t68 + 0.2e1 * t1 * t17 + 0.2e1 * t7 * t9 + 0.2e1 * t2 * t18 + 0.2e1 * t6 * t8 + ((mrSges(4,2) * t70 - t30 * t10 + t32 * t11 + t12 * t67 + t40 * t33) * qJD(3) + t53) * t33 + (t5 * t67 - t30 * t3 + t32 * t4 + (-t30 * t11 - t32 * t10 - t33 * (-t57 - t59)) * qJD(4) + (mrSges(4,1) * t70 + (Ifges(5,5) * t32 - t40) * t31 + (t24 ^ 2 * t68 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3)) * t33) * qJD(3)) * t31; -t33 * t5 + (m(5) * (t1 * t32 - t2 * t30 - t6 * t45 - t7 * t47) - t17 * t47 + t32 * t9 - t18 * t45 - t30 * t8) * t31 + (t31 * t12 + m(5) * (t7 * t54 - t6 * t55 + (t31 ^ 2 - t33 ^ 2) * t24) + t17 * t54 - t18 * t55) * qJD(3); (-0.1e1 + t52) * t31 * t48 * t68; -pkin(3) * t5 + (-t26 / 0.2e1 + (Ifges(4,5) + (-m(5) * pkin(3) + t69) * t24) * qJD(3)) * t33 + (-t2 * mrSges(5,3) + t4 / 0.2e1 + t48 * t66 + (-t10 / 0.2e1 + t56 / 0.2e1 - t7 * mrSges(5,3)) * qJD(4) + (m(5) * (-t2 - t50) - qJD(4) * t17 - t8) * pkin(6)) * t30 + (qJD(4) * t11 / 0.2e1 + t3 / 0.2e1 + t22 * t48 / 0.2e1 + t39 * mrSges(5,3) + (m(5) * t39 - qJD(4) * t18 + t9) * pkin(6)) * t32 + (t15 * t65 + t24 * t14 + t32 * t16 / 0.2e1 + (t22 * t65 + t32 * t66) * qJD(4) + (-Ifges(4,6) + t59 / 0.2e1 + t57 / 0.2e1 + t24 * mrSges(4,2)) * qJD(3)) * t31; (m(5) * (t52 * t63 - t64) + t69 * t31) * qJD(3) + (-t14 + (t52 * mrSges(5,3) - mrSges(4,2)) * qJD(3)) * t33; -0.2e1 * pkin(3) * t14 + t32 * t15 + t30 * t16 + (-t21 * t30 + t22 * t32) * qJD(4); t2 * mrSges(5,1) - t1 * mrSges(5,2) - Ifges(5,5) * t43 - t34 * Ifges(5,6) - t53; -t5; t26 + (t20 * pkin(6) - t58) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t23(1), t23(2), t23(4), t23(7); t23(2), t23(3), t23(5), t23(8); t23(4), t23(5), t23(6), t23(9); t23(7), t23(8), t23(9), t23(10);];
Mq = res;

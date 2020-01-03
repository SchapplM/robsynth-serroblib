% Calculate time derivative of joint inertia matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:09
% EndTime: 2019-12-31 16:56:10
% DurationCPUTime: 0.55s
% Computational Cost: add. (349->127), mult. (823->207), div. (0->0), fcn. (525->4), ass. (0->64)
t64 = 2 * qJD(2);
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t20 = -t30 * mrSges(5,1) + t28 * mrSges(5,2);
t67 = -mrSges(4,1) + t20;
t51 = t28 ^ 2 + t30 ^ 2;
t66 = 2 * m(5);
t32 = (-pkin(1) - pkin(5));
t65 = -2 * t32;
t63 = -t28 / 0.2e1;
t29 = sin(qJ(3));
t62 = t29 * pkin(3);
t31 = cos(qJ(3));
t61 = t31 * pkin(6);
t60 = Ifges(5,4) * t28;
t59 = Ifges(5,4) * t30;
t58 = Ifges(5,5) * t28;
t57 = Ifges(5,6) * t28;
t56 = Ifges(5,6) * t29;
t55 = Ifges(5,6) * t30;
t54 = t29 * t32;
t53 = t31 * mrSges(5,3);
t48 = qJD(3) * t29;
t42 = t28 * t48;
t47 = qJD(3) * t31;
t52 = Ifges(5,6) * t42 + Ifges(5,3) * t47;
t19 = qJ(2) - t61 + t62;
t8 = t30 * t19 - t28 * t54;
t50 = qJD(4) * t8;
t9 = t28 * t19 + t30 * t54;
t49 = qJD(4) * t9;
t46 = qJD(4) * t28;
t45 = qJD(4) * t30;
t44 = qJD(4) * t31;
t43 = t30 * t44;
t41 = t28 * t47;
t40 = t30 * t47;
t39 = -Ifges(5,5) * t30 + (2 * Ifges(4,4));
t13 = qJD(2) + (pkin(3) * t31 + pkin(6) * t29) * qJD(3);
t1 = t28 * t13 + t32 * t40 + t50;
t38 = t1 - t50;
t37 = mrSges(5,1) * t28 + mrSges(5,2) * t30;
t36 = Ifges(5,1) * t30 - t60;
t22 = Ifges(5,1) * t28 + t59;
t35 = -Ifges(5,2) * t28 + t59;
t21 = Ifges(5,2) * t30 + t60;
t34 = t28 * t44 + t30 * t48;
t33 = t42 - t43;
t25 = Ifges(5,5) * t45;
t18 = t29 * mrSges(5,1) - t30 * t53;
t17 = -t29 * mrSges(5,2) - t28 * t53;
t16 = t36 * qJD(4);
t15 = t35 * qJD(4);
t14 = t37 * qJD(4);
t12 = t37 * t31;
t11 = Ifges(5,5) * t29 + t36 * t31;
t10 = t35 * t31 + t56;
t7 = -mrSges(5,2) * t47 + t33 * mrSges(5,3);
t6 = mrSges(5,1) * t47 + t34 * mrSges(5,3);
t5 = -t33 * mrSges(5,1) - t34 * mrSges(5,2);
t4 = -t22 * t44 + (Ifges(5,5) * t31 - t36 * t29) * qJD(3);
t3 = -t21 * t44 + (Ifges(5,6) * t31 - t35 * t29) * qJD(3);
t2 = t30 * t13 - t32 * t41 - t49;
t23 = [(t9 * t1 + t8 * t2) * t66 + 0.2e1 * t1 * t17 + 0.2e1 * t9 * t7 + 0.2e1 * t2 * t18 + 0.2e1 * t8 * t6 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t64 + (mrSges(4,1) * t64 + (-0.2e1 * qJ(2) * mrSges(4,2) + t28 * t10 - t30 * t11 + 0.2e1 * t32 * t12 + t39 * t29) * qJD(3) + t52) * t29 + (mrSges(4,2) * t64 - t28 * t3 + t30 * t4 + t5 * t65 + (t29 * (-t55 - t58) - t28 * t11 - t30 * t10) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) + (-t39 - t57) * t31 + (-(2 * m(5) * t32 ^ 2) - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + Ifges(5,3)) * t29) * qJD(3)) * t31; (-t5 + (m(5) * (-t28 * t8 + t30 * t9) + t30 * t17 - t28 * t18) * qJD(3)) * t31 + (m(5) * (t1 * t30 - t2 * t28 - t8 * t45 - t9 * t46 + t47 * t65) - t17 * t46 + t30 * t7 + qJD(3) * t12 - t18 * t45 - t28 * t6) * t29; (-0.1e1 + t51) * t29 * t47 * t66; -pkin(3) * t5 + (t25 / 0.2e1 + (-Ifges(4,5) + (-m(5) * pkin(3) + t67) * t32) * qJD(3)) * t29 + (t4 / 0.2e1 - t2 * mrSges(5,3) + t21 * t48 / 0.2e1 + (-t56 / 0.2e1 - t10 / 0.2e1 - t9 * mrSges(5,3)) * qJD(4) + (m(5) * (-t2 - t49) - qJD(4) * t17 - t6) * pkin(6)) * t28 + (t3 / 0.2e1 + qJD(4) * t11 / 0.2e1 - t22 * t48 / 0.2e1 + t38 * mrSges(5,3) + (m(5) * t38 - qJD(4) * t18 + t7) * pkin(6)) * t30 + (-t32 * t14 + t30 * t16 / 0.2e1 + t15 * t63 + (-t30 * t21 / 0.2e1 + t22 * t63) * qJD(4) + (t58 / 0.2e1 + t55 / 0.2e1 - Ifges(4,6) - t32 * mrSges(4,2)) * qJD(3)) * t31; -t31 * t14 + (-t31 * mrSges(4,2) + m(5) * (t51 * t61 - t62) + t51 * t53 + t67 * t29) * qJD(3); -0.2e1 * pkin(3) * t14 + t30 * t15 + t28 * t16 + (-t21 * t28 + t22 * t30) * qJD(4); t2 * mrSges(5,1) - t1 * mrSges(5,2) - Ifges(5,5) * t34 - Ifges(5,6) * t43 + t52; (t29 * t46 - t40) * mrSges(5,2) + (-t29 * t45 - t41) * mrSges(5,1); t25 + (pkin(6) * t20 - t57) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t23(1), t23(2), t23(4), t23(7); t23(2), t23(3), t23(5), t23(8); t23(4), t23(5), t23(6), t23(9); t23(7), t23(8), t23(9), t23(10);];
Mq = res;

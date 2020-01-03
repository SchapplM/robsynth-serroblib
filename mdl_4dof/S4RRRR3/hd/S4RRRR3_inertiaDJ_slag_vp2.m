% Calculate time derivative of joint inertia matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:19
% DurationCPUTime: 0.66s
% Computational Cost: add. (1185->122), mult. (2736->214), div. (0->0), fcn. (2328->6), ass. (0->56)
t50 = sin(qJ(2));
t66 = -pkin(6) - pkin(5);
t44 = t66 * t50;
t53 = cos(qJ(2));
t45 = t66 * t53;
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t27 = t52 * t44 + t45 * t49;
t28 = t49 * t44 - t52 * t45;
t69 = qJD(2) + qJD(3);
t68 = 2 * m(5);
t47 = -pkin(2) * t53 - pkin(1);
t67 = 0.2e1 * t47;
t46 = pkin(2) * t52 + pkin(3);
t51 = cos(qJ(4));
t60 = qJD(4) * t51;
t48 = sin(qJ(4));
t61 = qJD(4) * t48;
t63 = t48 * t49;
t19 = t46 * t60 + (-t49 * t61 + (t51 * t52 - t63) * qJD(3)) * pkin(2);
t65 = t19 * mrSges(5,2);
t62 = t49 * t51;
t59 = 0.2e1 * t53;
t58 = qJD(2) * t66;
t20 = -t46 * t61 + (-t49 * t60 + (-t48 * t52 - t62) * qJD(3)) * pkin(2);
t18 = t20 * mrSges(5,1);
t57 = t18 - t65;
t38 = t49 * t53 + t52 * t50;
t15 = -pkin(7) * t38 + t27;
t37 = -t49 * t50 + t52 * t53;
t16 = pkin(7) * t37 + t28;
t10 = t15 * t51 - t16 * t48;
t42 = t50 * t58;
t43 = t53 * t58;
t13 = t27 * qJD(3) + t52 * t42 + t49 * t43;
t26 = t69 * t38;
t8 = -pkin(7) * t26 + t13;
t14 = -t28 * qJD(3) - t42 * t49 + t52 * t43;
t25 = t69 * t37;
t9 = -pkin(7) * t25 + t14;
t2 = t10 * qJD(4) + t48 * t9 + t51 * t8;
t11 = t15 * t48 + t16 * t51;
t3 = -t11 * qJD(4) - t48 * t8 + t51 * t9;
t21 = t37 * t51 - t38 * t48;
t6 = t21 * qJD(4) + t25 * t51 - t26 * t48;
t22 = t37 * t48 + t38 * t51;
t7 = -t22 * qJD(4) - t25 * t48 - t26 * t51;
t56 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t6 + Ifges(5,6) * t7;
t55 = (-mrSges(4,1) * t49 - mrSges(4,2) * t52) * qJD(3) * pkin(2);
t54 = t14 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,5) * t25 - Ifges(4,6) * t26 + t56;
t35 = (-t48 * mrSges(5,1) - t51 * mrSges(5,2)) * qJD(4) * pkin(3);
t31 = pkin(2) * t62 + t46 * t48;
t30 = -pkin(2) * t63 + t46 * t51;
t29 = -t37 * pkin(3) + t47;
t17 = pkin(2) * qJD(2) * t50 + pkin(3) * t26;
t1 = [(t26 * mrSges(4,1) + t25 * mrSges(4,2)) * t67 + 0.2e1 * t29 * (-t7 * mrSges(5,1) + t6 * mrSges(5,2)) - 0.2e1 * t37 * Ifges(4,2) * t26 + 0.2e1 * t25 * t38 * Ifges(4,1) + 0.2e1 * t21 * Ifges(5,2) * t7 + 0.2e1 * t17 * (-t21 * mrSges(5,1) + t22 * mrSges(5,2)) + 0.2e1 * t6 * t22 * Ifges(5,1) + 0.2e1 * m(4) * (t13 * t28 + t14 * t27) + (t10 * t3 + t11 * t2 + t17 * t29) * t68 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t53) * t59 + (-0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t67 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t37 + mrSges(4,2) * t38) - 0.2e1 * Ifges(3,4) * t50 + (Ifges(3,1) - Ifges(3,2)) * t59) * t50) * qJD(2) + 0.2e1 * (t21 * t6 + t7 * t22) * Ifges(5,4) + 0.2e1 * (t25 * t37 - t26 * t38) * Ifges(4,4) + 0.2e1 * (-t10 * t6 + t11 * t7 + t2 * t21 - t22 * t3) * mrSges(5,3) + 0.2e1 * (t13 * t37 - t14 * t38 - t25 * t27 - t26 * t28) * mrSges(4,3); m(5) * (t10 * t20 + t11 * t19 + t2 * t31 + t3 * t30) + (Ifges(3,5) * t53 - Ifges(3,6) * t50 + (-mrSges(3,1) * t53 + mrSges(3,2) * t50) * pkin(5)) * qJD(2) + (t19 * t21 - t20 * t22 - t30 * t6 + t31 * t7) * mrSges(5,3) + (m(4) * (t13 * t49 + t14 * t52 + (-t27 * t49 + t28 * t52) * qJD(3)) + (-t52 * t25 - t49 * t26 + (t37 * t52 + t38 * t49) * qJD(3)) * mrSges(4,3)) * pkin(2) + t54; (t19 * t31 + t20 * t30) * t68 - 0.2e1 * t65 + 0.2e1 * t18 + 0.2e1 * t55; (m(5) * (t2 * t48 + t3 * t51 + (-t10 * t48 + t11 * t51) * qJD(4)) + (t48 * t7 - t51 * t6 + (t21 * t51 + t22 * t48) * qJD(4)) * mrSges(5,3)) * pkin(3) + t54; t55 + (m(5) * (t19 * t48 + t20 * t51 - t30 * t61 + t31 * t60) - mrSges(5,2) * t60 - mrSges(5,1) * t61) * pkin(3) + t57; 0.2e1 * t35; t56; t57; t35; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S5RPRPR4
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
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:22:26
% DurationCPUTime: 0.46s
% Computational Cost: add. (1093->116), mult. (2295->194), div. (0->0), fcn. (2063->8), ass. (0->57)
t52 = cos(qJ(3));
t58 = -cos(pkin(8)) * pkin(1) - pkin(2);
t43 = -t52 * pkin(3) + t58;
t68 = 0.2e1 * t43;
t67 = m(5) * pkin(3);
t46 = sin(pkin(9));
t66 = pkin(3) * t46;
t47 = cos(pkin(9));
t50 = sin(qJ(3));
t41 = -t46 * t50 + t47 * t52;
t42 = t46 * t52 + t47 * t50;
t49 = sin(qJ(5));
t51 = cos(qJ(5));
t21 = t49 * t41 + t51 * t42;
t37 = t42 * qJD(3);
t38 = t41 * qJD(3);
t11 = -t21 * qJD(5) - t51 * t37 - t49 * t38;
t20 = t51 * t41 - t49 * t42;
t65 = t20 * t11;
t10 = t20 * qJD(5) - t49 * t37 + t51 * t38;
t64 = t21 * t10;
t63 = t41 * t37;
t62 = t42 * t38;
t44 = sin(pkin(8)) * pkin(1) + pkin(6);
t61 = qJ(4) + t44;
t56 = qJD(3) * t61;
t27 = t52 * qJD(4) - t50 * t56;
t28 = -t50 * qJD(4) - t52 * t56;
t15 = t47 * t27 + t46 * t28;
t39 = t61 * t50;
t40 = t61 * t52;
t19 = -t46 * t39 + t47 * t40;
t60 = 0.2e1 * t52;
t59 = qJD(3) * t50 * pkin(3);
t4 = -t11 * mrSges(6,1) + t10 * mrSges(6,2);
t57 = t37 * mrSges(5,1) + t38 * mrSges(5,2);
t14 = -t46 * t27 + t47 * t28;
t18 = -t47 * t39 - t46 * t40;
t12 = -t38 * pkin(7) + t14;
t13 = -t37 * pkin(7) + t15;
t16 = -t42 * pkin(7) + t18;
t17 = t41 * pkin(7) + t19;
t5 = t51 * t16 - t49 * t17;
t2 = t5 * qJD(5) + t49 * t12 + t51 * t13;
t6 = t49 * t16 + t51 * t17;
t3 = -t6 * qJD(5) + t51 * t12 - t49 * t13;
t55 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t11;
t45 = t47 * pkin(3) + pkin(4);
t35 = t51 * t45 - t49 * t66;
t29 = t35 * qJD(5);
t36 = t49 * t45 + t51 * t66;
t30 = t36 * qJD(5);
t54 = -t30 * mrSges(6,1) - t29 * mrSges(6,2);
t53 = -t4 - t57;
t26 = t37 * pkin(4) + t59;
t25 = -t41 * pkin(4) + t43;
t1 = [-0.2e1 * Ifges(5,2) * t63 + 0.2e1 * Ifges(5,1) * t62 + t57 * t68 + 0.2e1 * Ifges(6,1) * t64 + 0.2e1 * Ifges(6,2) * t65 + 0.2e1 * t25 * t4 + 0.2e1 * t26 * (-t20 * mrSges(6,1) + t21 * mrSges(6,2)) + 0.2e1 * m(5) * (t18 * t14 + t19 * t15) + 0.2e1 * m(6) * (t6 * t2 + t25 * t26 + t5 * t3) + ((t58 * mrSges(4,2) + Ifges(4,4) * t52) * t60 + (0.2e1 * pkin(3) * (-t41 * mrSges(5,1) + t42 * mrSges(5,2)) + 0.2e1 * t58 * mrSges(4,1) + t67 * t68 - 0.2e1 * Ifges(4,4) * t50 + (Ifges(4,1) - Ifges(4,2)) * t60) * t50) * qJD(3) + 0.2e1 * (t10 * t20 + t11 * t21) * Ifges(6,4) + 0.2e1 * (-t37 * t42 + t41 * t38) * Ifges(5,4) + 0.2e1 * (-t5 * t10 + t6 * t11 + t2 * t20 - t3 * t21) * mrSges(6,3) + 0.2e1 * (-t14 * t42 + t15 * t41 - t18 * t38 - t19 * t37) * mrSges(5,3); m(5) * (t14 * t41 + t15 * t42 - t18 * t37 + t19 * t38) + m(6) * (t6 * t10 + t5 * t11 + t2 * t21 + t3 * t20); 0.2e1 * m(5) * (t62 - t63) + 0.2e1 * m(6) * (t64 + t65); m(6) * (t36 * t2 + t29 * t6 + t35 * t3 - t30 * t5) - t15 * mrSges(5,2) + t14 * mrSges(5,1) + Ifges(5,5) * t38 - Ifges(5,6) * t37 + (Ifges(4,5) * t52 - Ifges(4,6) * t50 + (-mrSges(4,1) * t52 + mrSges(4,2) * t50) * t44) * qJD(3) + (m(5) * (t14 * t47 + t15 * t46) + (-t46 * t37 - t47 * t38) * mrSges(5,3)) * pkin(3) + (-t35 * t10 + t36 * t11 + t29 * t20 + t30 * t21) * mrSges(6,3) + t55; (-t50 * mrSges(4,1) - t52 * mrSges(4,2)) * qJD(3) + m(6) * (t36 * t10 + t35 * t11 - t30 * t20 + t29 * t21) + (-t37 * t47 + t38 * t46) * t67 + t53; 0.2e1 * m(6) * (t36 * t29 - t35 * t30) + 0.2e1 * t54; m(5) * t59 + m(6) * t26 - t53; 0; 0; 0; t55; -t4; t54; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

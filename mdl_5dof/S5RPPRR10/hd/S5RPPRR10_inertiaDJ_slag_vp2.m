% Calculate time derivative of joint inertia matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:00
% EndTime: 2019-12-31 18:04:01
% DurationCPUTime: 0.58s
% Computational Cost: add. (847->114), mult. (1936->194), div. (0->0), fcn. (1776->6), ass. (0->50)
t49 = sin(pkin(8));
t66 = -pkin(6) + qJ(2);
t37 = t66 * t49;
t50 = cos(pkin(8));
t38 = t66 * t50;
t52 = sin(qJ(4));
t54 = cos(qJ(4));
t22 = t54 * t37 - t38 * t52;
t70 = qJD(2) * (t49 ^ 2 + t50 ^ 2);
t23 = t52 * t37 + t54 * t38;
t62 = t49 * qJD(3);
t69 = qJD(4) + qJD(5);
t68 = 2 * m(6);
t65 = qJ(2) * t70;
t64 = qJD(2) * t52;
t63 = qJD(2) * t54;
t60 = -t50 * pkin(2) - t49 * qJ(3) - pkin(1);
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t35 = -t51 * t52 + t53 * t54;
t20 = t69 * t35;
t36 = t51 * t54 + t53 * t52;
t21 = t69 * t36;
t59 = -t21 * mrSges(6,1) - t20 * mrSges(6,2);
t27 = t50 * pkin(3) - t60;
t33 = -t49 * t52 - t50 * t54;
t34 = t49 * t54 - t50 * t52;
t16 = t33 * t53 - t34 * t51;
t25 = t33 * qJD(4);
t26 = t34 * qJD(4);
t8 = t16 * qJD(5) + t25 * t53 - t26 * t51;
t17 = t33 * t51 + t34 * t53;
t9 = -t17 * qJD(5) - t25 * t51 - t26 * t53;
t57 = -t9 * mrSges(6,1) + t8 * mrSges(6,2);
t12 = qJD(4) * t22 + t49 * t64 + t50 * t63;
t10 = -pkin(7) * t26 + t12;
t13 = -qJD(4) * t23 + t49 * t63 - t50 * t64;
t11 = -t25 * pkin(7) + t13;
t14 = -pkin(7) * t34 + t22;
t15 = pkin(7) * t33 + t23;
t4 = t14 * t53 - t15 * t51;
t2 = t4 * qJD(5) + t10 * t53 + t11 * t51;
t5 = t14 * t51 + t15 * t53;
t3 = -t5 * qJD(5) - t10 * t51 + t11 * t53;
t56 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t8 + Ifges(6,6) * t9;
t55 = t26 * mrSges(5,1) + t25 * mrSges(5,2);
t29 = (-mrSges(6,1) * t51 - mrSges(6,2) * t53) * qJD(5) * pkin(4);
t24 = pkin(4) * t26 + t62;
t19 = -pkin(4) * t33 + t27;
t1 = [0.2e1 * t27 * t55 - 0.2e1 * t33 * Ifges(5,2) * t26 + 0.2e1 * t34 * t25 * Ifges(5,1) + 0.2e1 * t16 * Ifges(6,2) * t9 + 0.2e1 * t8 * t17 * Ifges(6,1) + 0.2e1 * t19 * t57 + 0.2e1 * t24 * (-t16 * mrSges(6,1) + t17 * mrSges(6,2)) + 0.2e1 * (mrSges(4,1) * t50 - mrSges(5,1) * t33 + mrSges(5,2) * t34 + mrSges(4,3) * t49) * t62 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * t70 + 0.2e1 * (t16 * t8 + t17 * t9) * Ifges(6,4) + 0.2e1 * (t25 * t33 - t26 * t34) * Ifges(5,4) + 0.2e1 * (t16 * t2 - t17 * t3 - t4 * t8 + t5 * t9) * mrSges(6,3) + 0.2e1 * (t12 * t33 - t13 * t34 - t22 * t25 - t23 * t26) * mrSges(5,3) + 0.2e1 * m(3) * t65 + 0.2e1 * m(4) * (-t60 * t62 + t65) + 0.2e1 * m(5) * (t12 * t23 + t13 * t22 + t27 * t62) + (t19 * t24 + t2 * t5 + t3 * t4) * t68; -m(6) * t24 + (-m(4) - m(5)) * t62 - t55 - t57; 0; m(4) * t49 * qJD(2) + m(6) * (t2 * t36 + t20 * t5 - t21 * t4 + t3 * t35) + m(5) * (t52 * t12 + t13 * t54 + (-t22 * t52 + t23 * t54) * qJD(4)) + (t16 * t20 + t17 * t21 - t35 * t8 + t36 * t9) * mrSges(6,3) + (-t54 * t25 - t52 * t26 + (t33 * t54 + t34 * t52) * qJD(4)) * mrSges(5,3); 0; (t20 * t36 - t21 * t35) * t68; t13 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,5) * t25 - Ifges(5,6) * t26 + (m(6) * (t2 * t51 + t3 * t53 + (-t4 * t51 + t5 * t53) * qJD(5)) + (t51 * t9 - t53 * t8 + (t16 * t53 + t17 * t51) * qJD(5)) * mrSges(6,3)) * pkin(4) + t56; 0; (-mrSges(5,1) * t52 - mrSges(5,2) * t54) * qJD(4) + m(6) * (t20 * t51 - t21 * t53 + (-t35 * t51 + t36 * t53) * qJD(5)) * pkin(4) + t59; 0.2e1 * t29; t56; 0; t59; t29; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:08:40
% DurationCPUTime: 0.42s
% Computational Cost: add. (431->88), mult. (970->133), div. (0->0), fcn. (736->6), ass. (0->46)
t49 = (mrSges(6,2) + mrSges(5,3));
t57 = 2 * t49;
t56 = m(5) + m(6);
t26 = sin(pkin(8));
t52 = t26 * pkin(3);
t22 = qJ(5) + t52;
t55 = m(6) * t22 + mrSges(6,3);
t29 = cos(qJ(3));
t40 = -cos(pkin(7)) * pkin(1) - pkin(2);
t21 = -t29 * pkin(3) + t40;
t54 = 0.2e1 * t21;
t53 = m(5) * pkin(3);
t47 = cos(pkin(8));
t28 = sin(qJ(3));
t50 = t26 * t28;
t33 = t29 * t47 - t50;
t17 = t33 * qJD(3);
t51 = t17 * mrSges(6,2);
t24 = sin(pkin(7)) * pkin(1) + pkin(6);
t48 = qJ(4) + t24;
t46 = qJD(3) * t28;
t45 = qJD(3) * t29;
t37 = t47 * t28;
t20 = t26 * t29 + t37;
t44 = t20 * qJD(5);
t43 = 0.2e1 * t29;
t42 = pkin(3) * t46;
t36 = qJD(3) * t48;
t12 = t29 * qJD(4) - t28 * t36;
t30 = -t28 * qJD(4) - t29 * t36;
t5 = t26 * t12 - t30 * t47;
t6 = t12 * t47 + t26 * t30;
t18 = t48 * t29;
t8 = t26 * t18 + t37 * t48;
t9 = t18 * t47 - t48 * t50;
t41 = t8 * t5 + t9 * t6;
t39 = t47 * pkin(3);
t35 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t16 = t20 * qJD(3);
t14 = t16 * mrSges(6,1);
t15 = t17 * mrSges(5,2);
t32 = -t16 * mrSges(5,1) + t17 * mrSges(6,3) - t14 - t15;
t25 = -t39 - pkin(4);
t7 = -pkin(4) * t33 - t20 * qJ(5) + t21;
t3 = t16 * pkin(4) - t17 * qJ(5) + t42 - t44;
t1 = [0.2e1 * t3 * (-mrSges(6,1) * t33 - t20 * mrSges(6,3)) + t15 * t54 + 0.2e1 * t7 * t14 + 0.2e1 * m(5) * t41 + 0.2e1 * m(6) * (t7 * t3 + t41) + (-0.2e1 * t7 * mrSges(6,3) - t33 * t35 + t8 * t57 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t20) * t17 + (mrSges(5,1) * t54 + t20 * t35 - 0.2e1 * t49 * t9 - 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t33) * t16 + ((t40 * mrSges(4,2) + Ifges(4,4) * t29) * t43 + (t53 * t54 + 0.2e1 * t40 * mrSges(4,1) + 0.2e1 * pkin(3) * (-mrSges(5,1) * t33 + t20 * mrSges(5,2)) - 0.2e1 * Ifges(4,4) * t28 + (-Ifges(4,2) + Ifges(4,1)) * t43) * t28) * qJD(3) + (t5 * t20 + t33 * t6) * t57; t56 * (t8 * t16 + t9 * t17 + t6 * t20 - t33 * t5); 0.2e1 * t56 * (-t16 * t33 + t20 * t17); Ifges(4,5) * t45 - Ifges(4,6) * t46 + t25 * t51 + (m(6) * t9 + mrSges(6,2) * t33) * qJD(5) + (t26 * t53 - mrSges(5,2) + t55) * t6 + (m(6) * t25 - t47 * t53 - mrSges(5,1) - mrSges(6,1)) * t5 + (-mrSges(4,1) * t45 + mrSges(4,2) * t46) * t24 + (-mrSges(5,3) * t39 + Ifges(6,4) + Ifges(5,5)) * t17 + (-mrSges(6,2) * t22 - mrSges(5,3) * t52 - Ifges(5,6) + Ifges(6,6)) * t16; -mrSges(4,2) * t45 - mrSges(4,1) * t46 + (-t16 * t47 + t17 * t26) * t53 + m(6) * (t25 * t16 + t22 * t17 + t44) + t32; 0.2e1 * t55 * qJD(5); m(5) * t42 + m(6) * t3 - t32; 0; 0; 0; m(6) * t5 + t51; m(6) * t16; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

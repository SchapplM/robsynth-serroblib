% Calculate time derivative of joint inertia matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:23
% DurationCPUTime: 0.43s
% Computational Cost: add. (325->79), mult. (804->119), div. (0->0), fcn. (618->4), ass. (0->34)
t55 = mrSges(5,3) + mrSges(6,2);
t27 = sin(pkin(7));
t28 = cos(pkin(7));
t54 = qJD(2) * (t27 ^ 2 + t28 ^ 2);
t53 = -mrSges(5,1) - mrSges(6,1);
t52 = m(6) * qJD(5);
t30 = cos(qJ(4));
t43 = qJD(4) * t30;
t42 = t27 * qJD(3);
t51 = m(6) * qJ(5) + mrSges(6,3);
t50 = -m(6) * pkin(4) + t53;
t49 = -mrSges(5,2) + t51;
t29 = sin(qJ(4));
t48 = t28 * t29;
t47 = -pkin(6) + qJ(2);
t46 = qJ(2) * t54;
t45 = qJD(2) * t27;
t44 = qJD(4) * t29;
t38 = -t28 * pkin(2) - t27 * qJ(3) - pkin(1);
t16 = t47 * t27;
t17 = t47 * t28;
t34 = t30 * t16 - t29 * t17;
t14 = t27 * t29 + t28 * t30;
t4 = qJD(2) * t14 + qJD(4) * t34;
t9 = t29 * t16 + t30 * t17;
t5 = qJD(2) * t48 + qJD(4) * t9 - t30 * t45;
t37 = -t34 * t5 + t9 * t4;
t12 = t28 * pkin(3) - t38;
t15 = t27 * t30 - t48;
t11 = t27 * t43 - t28 * t44;
t10 = qJD(4) * t14;
t6 = t14 * pkin(4) - t15 * qJ(5) + t12;
t2 = t11 * pkin(4) + t10 * qJ(5) - t15 * qJD(5) + t42;
t1 = [0.2e1 * t12 * (t11 * mrSges(5,1) - t10 * mrSges(5,2)) + 0.2e1 * t2 * (t14 * mrSges(6,1) - t15 * mrSges(6,3)) + 0.2e1 * t6 * (t11 * mrSges(6,1) + t10 * mrSges(6,3)) + 0.2e1 * (t28 * mrSges(4,1) + t14 * mrSges(5,1) + t15 * mrSges(5,2) + t27 * mrSges(4,3)) * t42 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * t54 + 0.2e1 * m(3) * t46 + 0.2e1 * m(4) * (-t38 * t42 + t46) + 0.2e1 * m(5) * (t12 * t42 + t37) + 0.2e1 * m(6) * (t6 * t2 + t37) + 0.2e1 * (Ifges(6,3) + Ifges(5,2)) * t14 * t11 - 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t15 * t10 + 0.2e1 * (Ifges(6,5) - Ifges(5,4)) * (-t14 * t10 + t15 * t11) + 0.2e1 * t55 * (t10 * t34 - t9 * t11 - t4 * t14 + t5 * t15); -m(6) * t2 + t53 * t11 - (-mrSges(5,2) + mrSges(6,3)) * t10 + (-m(4) - m(5)) * t42; 0; m(4) * t45 + (m(5) + m(6)) * (t29 * t4 - t30 * t5 - t34 * t44 + t9 * t43) + t55 * (t30 * t10 - t29 * t11 + (-t14 * t30 + t15 * t29) * qJD(4)); 0; 0; t9 * t52 + (Ifges(6,6) - Ifges(5,6)) * t11 - (Ifges(6,4) + Ifges(5,5)) * t10 + (pkin(4) * t10 - qJ(5) * t11 - qJD(5) * t14) * mrSges(6,2) + t50 * t5 + t49 * t4; 0; t49 * t43 + (t50 * qJD(4) + t52) * t29; 0.2e1 * t51 * qJD(5); m(6) * t5 - t10 * mrSges(6,2); 0; m(6) * t44; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

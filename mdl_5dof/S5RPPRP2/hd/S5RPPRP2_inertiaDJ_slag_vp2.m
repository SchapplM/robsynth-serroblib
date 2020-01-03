% Calculate time derivative of joint inertia matrix for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:06
% DurationCPUTime: 0.26s
% Computational Cost: add. (376->61), mult. (852->95), div. (0->0), fcn. (692->6), ass. (0->29)
t47 = 2 * mrSges(6,2) + 2 * mrSges(5,3);
t46 = m(5) + m(6);
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t27 = sin(qJ(4));
t28 = cos(qJ(4));
t19 = t28 * t24 + t27 * t25;
t15 = t19 * qJD(4);
t13 = t15 * mrSges(6,1);
t18 = t27 * t24 - t28 * t25;
t14 = t18 * qJD(4);
t37 = t15 * mrSges(5,1) - t14 * mrSges(5,2);
t6 = pkin(4) * t15 + t14 * qJ(5) - t19 * qJD(5);
t45 = -m(6) * t6 - t14 * mrSges(6,3) - t13 - t37;
t44 = m(6) * qJ(5) + mrSges(6,3);
t43 = -0.2e1 * mrSges(6,3);
t20 = sin(pkin(7)) * pkin(1) + qJ(3);
t41 = pkin(6) + t20;
t16 = t41 * t24;
t17 = t41 * t25;
t32 = -t28 * t16 - t27 * t17;
t3 = -qJD(3) * t18 + qJD(4) * t32;
t9 = -t27 * t16 + t28 * t17;
t4 = qJD(3) * t19 + qJD(4) * t9;
t38 = t9 * t3 - t32 * t4;
t33 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t31 = -cos(pkin(7)) * pkin(1) - t25 * pkin(3) - pkin(2);
t7 = t18 * pkin(4) - t19 * qJ(5) + t31;
t1 = [0.2e1 * t31 * t37 + 0.2e1 * t7 * t13 - t15 * t9 * t47 + 0.2e1 * m(5) * t38 + 0.2e1 * m(6) * (t7 * t6 + t38) + (t15 * t33 + t4 * t47 + t6 * t43) * t19 + (0.2e1 * t6 * mrSges(6,1) - t3 * t47 + 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t15) * t18 - (t7 * t43 - t32 * t47 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t19 + t18 * t33) * t14 + 0.2e1 * (m(4) * t20 + mrSges(4,3)) * qJD(3) * (t24 ^ 2 + t25 ^ 2); t46 * (-t9 * t14 - t15 * t32 + t4 * t18 + t3 * t19); 0.2e1 * t46 * (-t19 * t14 + t18 * t15); -t45; 0; 0; m(6) * qJD(5) * t9 + (Ifges(6,6) - Ifges(5,6)) * t15 - (Ifges(6,4) + Ifges(5,5)) * t14 + (pkin(4) * t14 - qJ(5) * t15 - qJD(5) * t18) * mrSges(6,2) + (-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t4 + (-mrSges(5,2) + t44) * t3; t45; 0; 0.2e1 * t44 * qJD(5); m(6) * t4 - t14 * mrSges(6,2); m(6) * t15; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

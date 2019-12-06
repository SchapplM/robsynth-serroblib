% Calculate time derivative of joint inertia matrix for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:19
% EndTime: 2019-12-05 15:26:20
% DurationCPUTime: 0.16s
% Computational Cost: add. (102->42), mult. (290->70), div. (0->0), fcn. (206->6), ass. (0->24)
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t7 = t16 * mrSges(6,1) + t18 * mrSges(6,2);
t31 = mrSges(5,3) + t7;
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t29 = -t14 * t17 + t15 * t19;
t13 = t18 ^ 2;
t3 = t29 * qJD(2);
t5 = t14 * t19 + t15 * t17;
t1 = t5 * t3;
t25 = t16 ^ 2 + t13;
t24 = qJD(5) * t16;
t23 = qJD(5) * t18;
t22 = -t15 * pkin(2) - pkin(3);
t2 = t5 * qJD(2);
t21 = t25 * t2;
t10 = t14 * pkin(2) + qJ(4);
t20 = t5 * qJD(4) + t10 * t3;
t9 = -pkin(6) + t22;
t6 = -mrSges(6,1) * t23 + mrSges(6,2) * t24;
t4 = [0.2e1 * m(6) * (-t21 * t29 + t1) + 0.2e1 * (m(4) + m(5)) * (-t2 * t29 + t1); -t5 * t6 + (-t17 * mrSges(3,1) - t19 * mrSges(3,2)) * qJD(2) + (-mrSges(4,2) + t31) * t3 + (-t25 * mrSges(6,3) - mrSges(4,1) + mrSges(5,2)) * t2 + m(4) * (t14 * t3 - t15 * t2) * pkin(2) + m(5) * (t22 * t2 + t20) + m(6) * (t9 * t21 + t20); -0.2e1 * t13 * Ifges(6,4) * qJD(5) - 0.2e1 * t10 * t6 + 0.2e1 * ((m(5) + m(6)) * t10 + t31) * qJD(4) + 0.2e1 * (Ifges(6,4) * t16 + (-Ifges(6,1) + Ifges(6,2)) * t18) * t24; 0; 0; 0; 0.2e1 * (m(5) / 0.2e1 + m(6) * t25 / 0.2e1) * t2; 0; 0; 0; (-t16 * t2 + t23 * t29) * mrSges(6,2) + (t18 * t2 + t24 * t29) * mrSges(6,1); ((-mrSges(6,2) * t9 - Ifges(6,6)) * t18 + (-mrSges(6,1) * t9 - Ifges(6,5)) * t16) * qJD(5); t6; -t7 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;

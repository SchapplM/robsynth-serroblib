% Calculate time derivative of joint inertia matrix for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (173->40), mult. (370->70), div. (0->0), fcn. (286->4), ass. (0->25)
t18 = sin(pkin(6));
t19 = cos(pkin(6));
t21 = cos(qJ(4));
t32 = sin(qJ(4));
t23 = t32 * t18 - t21 * t19;
t26 = (t18 ^ 2 + t19 ^ 2) * qJD(3);
t6 = t23 * qJD(4);
t8 = -t21 * t18 - t32 * t19;
t35 = t8 * t6;
t7 = t8 * qJD(4);
t34 = t23 * t7;
t20 = -pkin(1) - qJ(3);
t33 = -pkin(5) + t20;
t29 = 2 * mrSges(5,3);
t28 = -t6 * mrSges(5,1) + t7 * mrSges(5,2);
t25 = t34 - t35;
t10 = t33 * t18;
t11 = t33 * t19;
t3 = -t32 * t10 + t21 * t11;
t24 = -t21 * t10 - t32 * t11;
t1 = t8 * qJD(3) + t3 * qJD(4);
t2 = t23 * qJD(3) + t24 * qJD(4);
t22 = t8 * t1 + t2 * t23 - t24 * t6 - t7 * t3;
t13 = t18 * pkin(3) + qJ(2);
t4 = [0.2e1 * t13 * t28 - 0.2e1 * Ifges(5,1) * t34 + 0.2e1 * Ifges(5,2) * t35 + 0.2e1 * m(5) * (t13 * qJD(2) - t1 * t24 + t3 * t2) + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t20 * t26) + t22 * t29 + 0.2e1 * (m(3) * qJ(2) + t18 * mrSges(4,1) - t8 * mrSges(5,1) + t19 * mrSges(4,2) - mrSges(5,2) * t23 + mrSges(3,3)) * qJD(2) + 0.2e1 * mrSges(4,3) * t26 + 0.2e1 * (-t23 * t6 + t7 * t8) * Ifges(5,4); -m(4) * t26 - m(5) * t22 + t25 * t29; -0.2e1 * m(5) * t25; (m(4) + m(5)) * qJD(2) + t28; 0; 0; t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t7 + Ifges(5,6) * t6; t7 * mrSges(5,1) + t6 * mrSges(5,2); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t4(1), t4(2), t4(4), t4(7); t4(2), t4(3), t4(5), t4(8); t4(4), t4(5), t4(6), t4(9); t4(7), t4(8), t4(9), t4(10);];
Mq = res;

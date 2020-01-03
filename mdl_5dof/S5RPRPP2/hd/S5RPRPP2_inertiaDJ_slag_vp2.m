% Calculate time derivative of joint inertia matrix for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:42
% DurationCPUTime: 0.25s
% Computational Cost: add. (173->78), mult. (389->100), div. (0->0), fcn. (198->4), ass. (0->34)
t36 = 2 * mrSges(6,3);
t20 = -cos(pkin(7)) * pkin(1) - pkin(2);
t35 = 0.2e1 * t20;
t13 = sin(qJ(3));
t25 = qJD(3) * t13;
t26 = t13 * qJ(4);
t16 = -t20 + t26;
t14 = cos(qJ(3));
t30 = t14 * pkin(3);
t6 = -t16 - t30;
t34 = -0.2e1 * t6;
t23 = t13 * qJD(4);
t19 = -pkin(3) * t25 + t23;
t24 = qJD(3) * t14;
t5 = qJ(4) * t24 + t19;
t33 = m(5) * t5;
t32 = m(5) + m(6);
t31 = pkin(3) + pkin(4);
t29 = -mrSges(5,2) + mrSges(6,3);
t9 = sin(pkin(7)) * pkin(1) + pkin(6);
t28 = qJ(5) - t9;
t27 = qJ(4) * t14;
t22 = 0.2e1 * t14;
t21 = -Ifges(5,5) + Ifges(4,4) - Ifges(6,4);
t8 = t28 * t14;
t18 = -t14 * mrSges(5,1) - t13 * mrSges(5,3);
t17 = (m(5) * t9 - t29) * t14;
t10 = mrSges(6,2) * t24;
t7 = t28 * t13;
t4 = t31 * t14 + t16;
t3 = -qJD(3) * t8 - t13 * qJD(5);
t2 = -t14 * qJD(5) + t28 * t25;
t1 = (-pkin(4) * t13 + t27) * qJD(3) + t19;
t11 = [t33 * t34 + 0.2e1 * m(6) * (t4 * t1 - t8 * t2 - t7 * t3) - 0.2e1 * t5 * t18 + 0.2e1 * t1 * (t14 * mrSges(6,1) + t13 * mrSges(6,2)) + 0.2e1 * t4 * t10 + (-t3 * t13 - t2 * t14) * t36 + ((mrSges(4,2) * t35 + mrSges(5,3) * t34 + t21 * t22 + t7 * t36) * t14 + (0.2e1 * t6 * mrSges(5,1) - 0.2e1 * t4 * mrSges(6,1) + mrSges(4,1) * t35 - 0.2e1 * t8 * mrSges(6,3) - 0.2e1 * t21 * t13 + (-Ifges(6,2) + Ifges(5,1) - Ifges(5,3) + Ifges(4,1) - Ifges(4,2) + Ifges(6,1)) * t22) * t13) * qJD(3); m(6) * (t2 * t13 - t3 * t14 + (-t13 * t7 - t14 * t8) * qJD(3)); 0; m(6) * (qJ(4) * t2 - qJD(4) * t8 - t3 * t31) - t3 * mrSges(6,1) + t2 * mrSges(6,2) + qJD(4) * t17 + ((-pkin(3) * mrSges(5,2) + mrSges(6,3) * t31 + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t14 + (t29 * qJ(4) - Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t13 + (t13 * mrSges(4,2) - t14 * mrSges(4,1) + m(5) * (-t26 - t30) + t18) * t9) * qJD(3); t10 + m(6) * t23 + t33 + (m(6) * t27 + (-mrSges(4,2) + mrSges(5,3)) * t14) * qJD(3) + (-m(6) * t31 - mrSges(4,1) - mrSges(5,1) - mrSges(6,1)) * t25; 0.2e1 * (t32 * qJ(4) + mrSges(6,2) + mrSges(5,3)) * qJD(4); m(6) * t3 + qJD(3) * t17; t32 * t25; 0; 0; m(6) * t1 - mrSges(6,1) * t25 + t10; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;

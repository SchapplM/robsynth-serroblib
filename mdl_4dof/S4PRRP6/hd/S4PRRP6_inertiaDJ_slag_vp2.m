% Calculate time derivative of joint inertia matrix for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:14
% DurationCPUTime: 0.19s
% Computational Cost: add. (87->48), mult. (278->73), div. (0->0), fcn. (140->4), ass. (0->23)
t13 = cos(qJ(3));
t10 = t13 ^ 2;
t11 = sin(qJ(3));
t26 = t11 ^ 2 + t10;
t6 = -t13 * mrSges(5,1) - t11 * mrSges(5,3);
t28 = -t13 * mrSges(4,1) + t11 * mrSges(4,2) + t6;
t14 = cos(qJ(2));
t23 = qJD(2) * t14;
t27 = t26 * pkin(5) * t23;
t25 = Ifges(4,4) - Ifges(5,5);
t12 = sin(qJ(2));
t24 = qJD(2) * t12;
t21 = (m(5) * pkin(5) + mrSges(5,2)) * t13;
t19 = t11 * mrSges(4,1) + t13 * mrSges(4,2);
t18 = t11 * mrSges(5,1) - t13 * mrSges(5,3);
t17 = -t13 * pkin(3) - t11 * qJ(4);
t16 = pkin(3) * t11 - qJ(4) * t13;
t15 = m(5) * t17 + t28;
t5 = -pkin(2) + t17;
t4 = t19 * qJD(3);
t3 = t18 * qJD(3);
t2 = t16 * qJD(3) - t11 * qJD(4);
t1 = [0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t26) * t12 * t23; (-mrSges(3,1) + t28) * t24 + m(4) * (-pkin(2) * t24 + t27) + m(5) * (t5 * t24 + t27) + (-t4 - t3 - m(5) * t2 + (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t26) * qJD(2)) * t14; -0.2e1 * pkin(2) * t4 + 0.2e1 * t5 * t3 + 0.2e1 * (m(5) * t5 + t6) * t2 + 0.2e1 * (t25 * t10 + (-t25 * t11 + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t13) * t11) * qJD(3); (-m(5) * t16 - t18 - t19) * t23 + (m(5) * qJD(4) * t13 + t15 * qJD(3)) * t12; qJD(4) * t21 + ((-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t13 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t11 + t15 * pkin(5)) * qJD(3); 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); (qJD(3) * t13 * t12 + t11 * t23) * m(5); qJD(3) * t21; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

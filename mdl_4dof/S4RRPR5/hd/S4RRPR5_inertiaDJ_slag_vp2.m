% Calculate time derivative of joint inertia matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (103->45), mult. (271->63), div. (0->0), fcn. (114->4), ass. (0->24)
t14 = cos(qJ(4));
t35 = t14 ^ 2;
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t25 = pkin(1) * qJD(2);
t12 = sin(qJ(4));
t26 = t12 ^ 2 + t35;
t34 = (-mrSges(3,2) * t15 + (-t26 * mrSges(5,3) - mrSges(3,1) + mrSges(4,2)) * t13) * t25;
t4 = t12 * mrSges(5,1) + t14 * mrSges(5,2);
t33 = mrSges(4,3) + t4;
t20 = mrSges(5,1) * t14 - mrSges(5,2) * t12;
t2 = t20 * qJD(4);
t31 = 0.2e1 * t2;
t7 = t15 * t25 + qJD(3);
t9 = t13 * pkin(1) + qJ(3);
t30 = t9 * t7;
t23 = t13 * t25;
t22 = -t15 * pkin(1) - pkin(2);
t19 = qJ(3) * t7 + qJD(3) * t9;
t18 = t26 * t23;
t17 = (-0.2e1 * t35 * Ifges(5,4) + (0.2e1 * Ifges(5,4) * t12 + 0.2e1 * (-Ifges(5,1) + Ifges(5,2)) * t14) * t12) * qJD(4);
t16 = -pkin(2) - pkin(6);
t8 = -pkin(6) + t22;
t1 = [t9 * t31 + 0.2e1 * t33 * t7 + 0.2e1 * t34 + 0.2e1 * m(5) * (t8 * t18 + t30) + 0.2e1 * m(4) * (t22 * t23 + t30) + t17; (qJ(3) + t9) * t2 + t34 + m(5) * (t16 * t18 + t19) + m(4) * (-pkin(2) * t23 + t19) + t17 + t33 * (qJD(3) + t7); qJ(3) * t31 + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t33) * qJD(3) + t17; (m(5) * t26 + m(4)) * t23; 0; 0; t20 * t23 + ((-mrSges(5,2) * t8 - Ifges(5,6)) * t14 + (-mrSges(5,1) * t8 - Ifges(5,5)) * t12) * qJD(4); ((-mrSges(5,2) * t16 - Ifges(5,6)) * t14 + (-mrSges(5,1) * t16 - Ifges(5,5)) * t12) * qJD(4); -t4 * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

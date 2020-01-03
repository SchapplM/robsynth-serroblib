% Calculate time derivative of joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:24
% DurationCPUTime: 0.35s
% Computational Cost: add. (276->75), mult. (675->120), div. (0->0), fcn. (495->4), ass. (0->33)
t35 = (mrSges(5,2) + mrSges(4,3));
t37 = 2 * t35;
t25 = cos(qJ(2));
t21 = -t25 * pkin(2) - pkin(1);
t36 = 0.2e1 * t21;
t34 = -qJ(3) - pkin(5);
t33 = 0.2e1 * t25;
t24 = sin(qJ(2));
t32 = qJD(2) * t24 * pkin(2);
t29 = qJD(2) * t34;
t11 = t25 * qJD(3) + t24 * t29;
t22 = sin(pkin(6));
t23 = cos(pkin(6));
t26 = -t24 * qJD(3) + t25 * t29;
t3 = t22 * t11 - t23 * t26;
t4 = t23 * t11 + t22 * t26;
t17 = t34 * t25;
t30 = t34 * t24;
t6 = -t22 * t17 - t23 * t30;
t7 = -t23 * t17 + t22 * t30;
t31 = t6 * t3 + t7 * t4;
t28 = -2 * Ifges(4,4) + 2 * Ifges(5,5);
t15 = t22 * t25 + t23 * t24;
t14 = t22 * t24 - t23 * t25;
t20 = -t23 * pkin(2) - pkin(3);
t18 = t22 * pkin(2) + qJ(4);
t13 = t14 * qJD(2);
t12 = t15 * qJD(2);
t10 = t13 * mrSges(4,2);
t9 = t12 * mrSges(5,1);
t5 = t14 * pkin(3) - t15 * qJ(4) + t21;
t2 = t12 * pkin(3) + t13 * qJ(4) - t15 * qJD(4) + t32;
t1 = [-t10 * t36 + 0.2e1 * t5 * t9 + 0.2e1 * t2 * (t14 * mrSges(5,1) - t15 * mrSges(5,3)) + 0.2e1 * m(4) * t31 + 0.2e1 * m(5) * (t5 * t2 + t31) - (-0.2e1 * t5 * mrSges(5,3) + t14 * t28 + t6 * t37 + 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t15) * t13 + (mrSges(4,1) * t36 + t15 * t28 - 0.2e1 * t35 * t7 + 0.2e1 * (Ifges(4,2) + Ifges(5,3)) * t14) * t12 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t25) * t33 + (m(4) * pkin(2) * t36 - 0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (t14 * mrSges(4,1) + t15 * mrSges(4,2)) - 0.2e1 * Ifges(3,4) * t24 + (Ifges(3,1) - Ifges(3,2)) * t33) * t24) * qJD(2) + (-t4 * t14 + t3 * t15) * t37; m(5) * (qJD(4) * t7 + t18 * t4 + t20 * t3) - t4 * mrSges(4,2) - t3 * mrSges(4,1) - t3 * mrSges(5,1) + t4 * mrSges(5,3) - (Ifges(5,4) + Ifges(4,5)) * t13 + (Ifges(5,6) - Ifges(4,6)) * t12 + (-qJD(4) * t14 - t18 * t12 - t20 * t13) * mrSges(5,2) + (Ifges(3,5) * t25 - Ifges(3,6) * t24 + (-mrSges(3,1) * t25 + mrSges(3,2) * t24) * pkin(5)) * qJD(2) + (m(4) * (t22 * t4 - t23 * t3) + (-t22 * t12 + t23 * t13) * mrSges(4,3)) * pkin(2); 0.2e1 * (m(5) * t18 + mrSges(5,3)) * qJD(4); m(4) * t32 + m(5) * t2 + t12 * mrSges(4,1) + t13 * mrSges(5,3) - t10 + t9; 0; 0; m(5) * t3 - t13 * mrSges(5,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

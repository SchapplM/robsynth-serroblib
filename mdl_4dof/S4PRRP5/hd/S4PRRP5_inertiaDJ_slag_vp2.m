% Calculate time derivative of joint inertia matrix for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:28:52
% DurationCPUTime: 0.24s
% Computational Cost: add. (106->57), mult. (318->92), div. (0->0), fcn. (171->4), ass. (0->28)
t35 = m(5) * pkin(3);
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t27 = t14 ^ 2 + t16 ^ 2;
t17 = cos(qJ(2));
t26 = qJD(2) * t17;
t9 = -t16 * pkin(3) - pkin(2);
t34 = m(5) * t9 - t16 * mrSges(5,1) + t14 * mrSges(5,2);
t25 = qJD(3) * t14;
t23 = qJD(3) * t16;
t4 = mrSges(5,1) * t25 + mrSges(5,2) * t23;
t33 = t25 * t35 + t4;
t32 = -2 * mrSges(5,3);
t30 = mrSges(4,2) + mrSges(5,2);
t29 = Ifges(4,4) + Ifges(5,4);
t28 = -qJ(4) - pkin(5);
t15 = sin(qJ(2));
t24 = qJD(3) * t15;
t22 = 0.2e1 * t14;
t20 = mrSges(5,1) + t35;
t19 = qJD(3) * t28;
t18 = -mrSges(4,1) - t20;
t8 = t28 * t16;
t6 = t28 * t14;
t5 = (mrSges(4,1) * t14 + mrSges(4,2) * t16) * qJD(3);
t3 = -t14 * qJD(4) + t16 * t19;
t2 = t16 * qJD(4) + t14 * t19;
t1 = [0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t27) * t15 * t26; m(5) * t8 * t14 * t24 + m(4) * t27 * pkin(5) * t26 + (m(5) * (-t14 * t3 + t16 * t2 - t23 * t6) + (-m(4) * pkin(2) - t16 * mrSges(4,1) + t14 * mrSges(4,2) - mrSges(3,1) + t34) * qJD(2)) * t15 + (-t5 + (-mrSges(3,2) + m(5) * (-t14 * t6 - t16 * t8) + (mrSges(4,3) + mrSges(5,3)) * t27) * qJD(2) - t33) * t17; 0.2e1 * m(5) * (-t8 * t2 + t6 * t3) - 0.2e1 * pkin(2) * t5 + 0.2e1 * t9 * t4 + (t3 * t32 + (0.2e1 * t34 * pkin(3) - t29 * t22 - t8 * t32) * qJD(3)) * t14 + (0.2e1 * t2 * mrSges(5,3) + (t6 * t32 + 0.2e1 * t29 * t16 + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,2)) * t22) * qJD(3)) * t16; (t30 * t14 + t18 * t16) * t24 + (t18 * t14 - t30 * t16) * t26; -t2 * mrSges(5,2) + t20 * t3 + ((mrSges(4,2) * pkin(5) - Ifges(4,6) - Ifges(5,6)) * t14 + (-mrSges(4,1) * pkin(5) - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t16) * qJD(3); 0; m(5) * qJD(2) * t15; t33; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

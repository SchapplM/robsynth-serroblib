% Calculate time derivative of joint inertia matrix for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:11
% EndTime: 2019-12-31 17:00:11
% DurationCPUTime: 0.24s
% Computational Cost: add. (117->67), mult. (294->92), div. (0->0), fcn. (132->2), ass. (0->26)
t28 = 2 * mrSges(5,1);
t27 = -2 * pkin(1);
t12 = cos(qJ(2));
t11 = sin(qJ(2));
t21 = t11 * qJ(3);
t14 = -t12 * pkin(2) - t21;
t6 = -pkin(1) + t14;
t25 = -0.2e1 * t6;
t24 = pkin(3) + pkin(5);
t23 = -mrSges(4,1) - mrSges(5,1);
t10 = -pkin(2) - qJ(4);
t22 = qJ(3) * t12;
t20 = qJD(2) * t11;
t19 = qJ(3) * qJD(3);
t18 = 0.2e1 * t12;
t17 = Ifges(4,6) + Ifges(3,4) - Ifges(5,6);
t8 = t24 * t12;
t16 = pkin(2) * t20 - t11 * qJD(3);
t15 = t12 * mrSges(4,2) - t11 * mrSges(4,3);
t13 = (m(4) * pkin(5) - t23) * t12;
t7 = t24 * t11;
t5 = qJD(2) * t8;
t4 = t24 * t20;
t3 = t10 * t12 - pkin(1) - t21;
t1 = -t12 * qJD(4) + (qJ(4) * t11 - t22) * qJD(2) + t16;
t2 = [0.2e1 * t1 * (-t11 * mrSges(5,2) - t12 * mrSges(5,3)) + 0.2e1 * m(5) * (t3 * t1 - t8 * t4 + t7 * t5) + (t5 * t11 - t4 * t12) * t28 + (((mrSges(3,2) * t27) - 0.2e1 * t3 * mrSges(5,2) + mrSges(4,3) * t25 + t17 * t18 + t7 * t28) * t12 + ((mrSges(3,1) * t27) - 0.2e1 * t8 * mrSges(5,1) + mrSges(4,2) * t25 + 0.2e1 * t3 * mrSges(5,3) - 0.2e1 * t17 * t11 + (Ifges(3,1) - Ifges(3,2) + Ifges(4,2) - Ifges(5,2) - Ifges(4,3) + Ifges(5,3)) * t18) * t11) * qJD(2) + 0.2e1 * (m(4) * t6 + t15) * (-qJD(2) * t22 + t16); -qJD(4) * t11 * mrSges(5,1) + m(5) * (-qJ(3) * t4 + qJD(3) * t8 - qJD(4) * t7 + t10 * t5) - t4 * mrSges(5,2) - t5 * mrSges(5,3) + qJD(3) * t13 + ((-pkin(2) * mrSges(4,1) + t10 * mrSges(5,1) - Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t12 + (t23 * qJ(3) + Ifges(5,4) + Ifges(4,5) - Ifges(3,6)) * t11 + (m(4) * t14 - t12 * mrSges(3,1) + t11 * mrSges(3,2) + t15) * pkin(5)) * qJD(2); 0.2e1 * m(4) * t19 + 0.2e1 * qJD(4) * mrSges(5,3) + 0.2e1 * m(5) * (-t10 * qJD(4) + t19) + 0.2e1 * (mrSges(4,3) + mrSges(5,2)) * qJD(3); m(5) * t5 + qJD(2) * t13; -m(5) * qJD(4); 0; -m(5) * t4 - mrSges(5,1) * t20; m(5) * qJD(3); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;

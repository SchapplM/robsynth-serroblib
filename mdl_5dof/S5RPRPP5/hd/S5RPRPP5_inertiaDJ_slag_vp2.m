% Calculate time derivative of joint inertia matrix for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:09
% DurationCPUTime: 0.27s
% Computational Cost: add. (178->81), mult. (371->101), div. (0->0), fcn. (169->2), ass. (0->32)
t34 = 2 * mrSges(6,3);
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t37 = t11 * mrSges(4,1) + t12 * mrSges(4,2);
t13 = -pkin(3) - pkin(4);
t26 = t12 * qJ(4);
t5 = t11 * t13 - qJ(2) + t26;
t36 = -0.2e1 * t5;
t16 = -t11 * pkin(3) + t26;
t35 = 0.2e1 * qJ(2) - 0.2e1 * t16;
t33 = m(5) + m(6);
t23 = qJD(3) * t12;
t32 = qJ(4) * t23 + qJD(4) * t11;
t29 = -mrSges(5,2) + mrSges(6,3);
t28 = mrSges(6,2) + mrSges(5,3);
t27 = qJ(4) * t11;
t14 = (-pkin(1) - pkin(6));
t25 = qJ(5) + t14;
t24 = qJD(3) * t11;
t22 = 0.2e1 * t12;
t21 = -Ifges(5,5) - Ifges(6,4) + Ifges(4,4);
t20 = qJD(3) * t25;
t19 = -t12 * qJD(4) + qJD(2);
t2 = -t12 * qJD(5) + t11 * t20;
t3 = t11 * qJD(5) + t12 * t20;
t18 = t11 * t3 - t12 * t2;
t17 = t11 * mrSges(5,1) - t12 * mrSges(5,3);
t15 = (m(5) * t14 + t29) * t11;
t7 = t25 * t12;
t6 = t25 * t11;
t1 = (t13 * t12 - t27) * qJD(3) - t19;
t4 = [0.2e1 * m(6) * (t5 * t1 - t7 * t2 + t6 * t3) + 0.2e1 * t1 * (-t11 * mrSges(6,1) + t12 * mrSges(6,2)) + t18 * t34 + 0.2e1 * (mrSges(3,3) + (m(3) + m(4)) * qJ(2) + t37) * qJD(2) + ((0.2e1 * qJ(2) * mrSges(4,1) + mrSges(5,1) * t35 + mrSges(6,1) * t36 - t21 * t22 + t6 * t34) * t12 + (-0.2e1 * qJ(2) * mrSges(4,2) + mrSges(6,2) * t36 + mrSges(5,3) * t35 - t7 * t34 + 0.2e1 * t21 * t11 + (-Ifges(4,1) - Ifges(5,1) - Ifges(6,1) + Ifges(4,2) + Ifges(6,2) + Ifges(5,3)) * t22) * t11) * qJD(3) + (m(5) * t35 + 0.2e1 * t17) * ((pkin(3) * t12 + t27) * qJD(3) + t19); m(6) * ((-t7 * t11 + t6 * t12) * qJD(3) + t18); 0; m(6) * (qJ(4) * t3 + qJD(4) * t6 + t13 * t2) - t2 * mrSges(6,1) + t3 * mrSges(6,2) + qJD(4) * t15 + ((t29 * qJ(4) - Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t12 + (pkin(3) * mrSges(5,2) + t13 * mrSges(6,3) - Ifges(5,4) - Ifges(4,5) + Ifges(6,5)) * t11 + (m(5) * t16 - t17 - t37) * t14) * qJD(3); (-mrSges(4,1) - mrSges(5,1) - mrSges(6,1)) * t24 + m(5) * (-pkin(3) * t24 + t32) + m(6) * (t13 * t24 + t32) + (-mrSges(4,2) + t28) * t23; 0.2e1 * (t33 * qJ(4) + t28) * qJD(4); m(6) * t2 + qJD(3) * t15; t33 * t24; 0; 0; m(6) * t1 + (-t12 * mrSges(6,1) - t11 * mrSges(6,2)) * qJD(3); 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;

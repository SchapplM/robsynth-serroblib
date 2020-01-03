% Calculate time derivative of joint inertia matrix for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:40:48
% DurationCPUTime: 0.24s
% Computational Cost: add. (142->76), mult. (358->98), div. (0->0), fcn. (167->2), ass. (0->31)
t32 = 2 * mrSges(6,3);
t11 = sin(qJ(3));
t22 = qJD(3) * t11;
t31 = -2 * pkin(2);
t12 = cos(qJ(3));
t23 = t11 * qJ(4);
t15 = -t12 * pkin(3) - t23;
t6 = -pkin(2) + t15;
t30 = -0.2e1 * t6;
t20 = t11 * qJD(4);
t17 = -pkin(3) * t22 + t20;
t21 = qJD(3) * t12;
t3 = qJ(4) * t21 + t17;
t29 = m(5) * t3;
t28 = m(5) + m(6);
t27 = pkin(3) + pkin(4);
t26 = -mrSges(5,2) + mrSges(6,3);
t25 = pkin(6) - qJ(5);
t24 = qJ(4) * t12;
t19 = 0.2e1 * t12;
t18 = Ifges(5,5) - Ifges(4,4) + Ifges(6,4);
t8 = t25 * t12;
t16 = -t12 * mrSges(5,1) - t11 * mrSges(5,3);
t14 = (m(5) * pkin(6) - t26) * t12;
t9 = mrSges(6,2) * t21;
t7 = t25 * t11;
t5 = t27 * t12 + pkin(2) + t23;
t4 = qJD(3) * t8 - t11 * qJD(5);
t2 = -t12 * qJD(5) - t25 * t22;
t1 = (-pkin(4) * t11 + t24) * qJD(3) + t17;
t10 = [0; m(6) * (t2 * t11 - t4 * t12 + (t11 * t7 + t12 * t8) * qJD(3)); t29 * t30 + 0.2e1 * m(6) * (t5 * t1 + t8 * t2 + t7 * t4) - 0.2e1 * t3 * t16 + 0.2e1 * t1 * (t12 * mrSges(6,1) + t11 * mrSges(6,2)) + 0.2e1 * t5 * t9 + (-t4 * t11 - t2 * t12) * t32 + (((mrSges(4,2) * t31) + mrSges(5,3) * t30 - 0.2e1 * t7 * mrSges(6,3) - t18 * t19) * t12 + ((mrSges(4,1) * t31) + 0.2e1 * t6 * mrSges(5,1) - 0.2e1 * t5 * mrSges(6,1) + t8 * t32 + 0.2e1 * t18 * t11 + (Ifges(4,1) + Ifges(5,1) + Ifges(6,1) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3)) * t19) * t11) * qJD(3); t9 + m(6) * t20 + t29 + (m(6) * t24 + (-mrSges(4,2) + mrSges(5,3)) * t12) * qJD(3) + (-m(6) * t27 - mrSges(4,1) - mrSges(5,1) - mrSges(6,1)) * t22; m(6) * (qJ(4) * t2 + qJD(4) * t8 - t27 * t4) - t4 * mrSges(6,1) + t2 * mrSges(6,2) + qJD(4) * t14 + ((-pkin(3) * mrSges(5,2) + mrSges(6,3) * t27 + Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t12 + (t26 * qJ(4) - Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t11 + (m(5) * t15 - t12 * mrSges(4,1) + t11 * mrSges(4,2) + t16) * pkin(6)) * qJD(3); 0.2e1 * (t28 * qJ(4) + mrSges(6,2) + mrSges(5,3)) * qJD(4); t28 * t22; m(6) * t4 + qJD(3) * t14; 0; 0; 0; m(6) * t1 - mrSges(6,1) * t22 + t9; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;

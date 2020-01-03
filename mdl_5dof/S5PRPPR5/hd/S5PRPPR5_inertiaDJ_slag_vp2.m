% Calculate time derivative of joint inertia matrix for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:04
% DurationCPUTime: 0.23s
% Computational Cost: add. (139->59), mult. (378->97), div. (0->0), fcn. (257->6), ass. (0->31)
t20 = cos(qJ(5));
t15 = t20 ^ 2;
t18 = sin(qJ(5));
t29 = t18 ^ 2 + t15;
t36 = -t29 * mrSges(6,3) + mrSges(5,2);
t17 = cos(pkin(8));
t19 = sin(qJ(2));
t28 = qJD(2) * t19;
t16 = sin(pkin(8));
t21 = cos(qJ(2));
t31 = t21 * t16;
t1 = qJD(2) * t31 - t17 * t28;
t3 = t19 * t16 + t21 * t17;
t35 = t3 * t1;
t34 = t16 * t3;
t33 = t17 * t1;
t10 = t20 * mrSges(6,1) - t18 * mrSges(6,2);
t30 = mrSges(5,1) + t10;
t22 = -pkin(2) - pkin(3);
t9 = t17 * qJ(3) + t16 * t22;
t27 = qJD(5) * t18;
t2 = t3 * qJD(2);
t26 = t29 * t2;
t4 = t19 * t17 - t31;
t25 = t29 * t4;
t23 = -mrSges(6,1) * t18 - mrSges(6,2) * t20;
t8 = -t16 * qJ(3) + t17 * t22;
t7 = t23 * qJD(5);
t6 = -pkin(6) + t9;
t5 = pkin(4) - t8;
t11 = [0.2e1 * m(5) * (t4 * t2 + t35) + 0.2e1 * m(6) * (t2 * t25 + t35); t3 * t7 + t30 * t1 + t36 * t2 + ((-mrSges(3,2) + mrSges(4,3)) * t21 + (-mrSges(3,1) - mrSges(4,1)) * t19) * qJD(2) + m(4) * (t19 * qJD(3) + (-pkin(2) * t19 + qJ(3) * t21) * qJD(2)) + m(5) * (-t8 * t1 + t9 * t2 + (t17 * t4 + t34) * qJD(3)) + m(6) * (t5 * t1 + t6 * t26 + (t17 * t25 + t34) * qJD(3)); 0.2e1 * qJD(5) * t15 * Ifges(6,4) + 0.2e1 * t5 * t7 + 0.2e1 * (-Ifges(6,4) * t18 + (Ifges(6,1) - Ifges(6,2)) * t20) * t27 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3) + (-m(5) * t8 + m(6) * t5 + t30) * t16 + (m(6) * t29 * t6 + m(5) * t9 + t36) * t17) * qJD(3); m(4) * t28 + m(5) * (t16 * t2 - t33) + m(6) * (t16 * t26 - t33); (-t7 + m(6) * (-0.1e1 + t29) * t16 * qJD(3)) * t17; 0; 0; 0; 0; 0; (-t2 * t20 + t4 * t27) * mrSges(6,2) + (-qJD(5) * t20 * t4 - t18 * t2) * mrSges(6,1); t23 * t17 * qJD(3) + ((-mrSges(6,1) * t6 - Ifges(6,5)) * t20 + (mrSges(6,2) * t6 + Ifges(6,6)) * t18) * qJD(5); -t10 * t16 * qJD(5); t7; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;

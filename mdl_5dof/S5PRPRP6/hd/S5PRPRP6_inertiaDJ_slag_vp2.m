% Calculate time derivative of joint inertia matrix for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:22
% EndTime: 2019-12-05 15:40:23
% DurationCPUTime: 0.25s
% Computational Cost: add. (123->63), mult. (348->89), div. (0->0), fcn. (165->4), ass. (0->30)
t17 = cos(qJ(4));
t14 = t17 ^ 2;
t15 = sin(qJ(4));
t29 = -t15 ^ 2 - t14;
t10 = t15 * mrSges(5,1) + t17 * mrSges(5,2);
t37 = mrSges(4,3) + t10;
t36 = m(6) * qJD(5) * t15;
t35 = m(6) * qJ(5) + mrSges(6,3);
t34 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t33 = -mrSges(5,2) + t35;
t19 = -pkin(2) - pkin(6);
t16 = sin(qJ(2));
t28 = qJD(2) * t16;
t32 = t29 * t19 * t28;
t31 = Ifges(5,4) - Ifges(6,5);
t18 = cos(qJ(2));
t30 = qJD(2) * t18 * qJ(3) + t16 * qJD(3);
t27 = qJD(4) * t15;
t26 = qJD(4) * t18;
t22 = -t15 * pkin(4) + t17 * qJ(5);
t6 = qJ(3) - t22;
t9 = t15 * mrSges(6,1) - t17 * mrSges(6,3);
t25 = m(6) * t6 + t9;
t23 = (m(6) * t19 - mrSges(6,2)) * t15;
t21 = (m(5) / 0.2e1 + m(6) / 0.2e1) * t28;
t20 = m(6) * t22 - t10 - t9;
t5 = (mrSges(5,1) * t17 - mrSges(5,2) * t15) * qJD(4);
t4 = (mrSges(6,1) * t17 + mrSges(6,3) * t15) * qJD(4);
t2 = -t17 * qJD(5) + qJD(3) + (pkin(4) * t17 + qJ(5) * t15) * qJD(4);
t1 = [0.4e1 * (0.1e1 + t29) * t18 * t21; (t4 + t5) * t16 + m(4) * t30 + m(5) * (t30 - t32) + m(6) * (t2 * t16 - t32) + ((-mrSges(3,2) + t25 + t37) * t18 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t29) * t16) * qJD(2); 0.2e1 * qJ(3) * t5 + 0.2e1 * t6 * t4 + 0.2e1 * t25 * t2 - 0.2e1 * t31 * qJD(4) * t14 + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t37) * qJD(3) + 0.2e1 * (t31 * t15 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,3)) * t17) * t27; m(4) * t28 - 0.2e1 * t29 * t21; 0; 0; -t18 * t36 + (t33 * t15 + t34 * t17) * t28 + (t34 * t15 - t33 * t17) * t26; qJD(5) * t23 + ((-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t17 + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t15 + t20 * t19) * qJD(4); t20 * qJD(4) + t36; 0.2e1 * t35 * qJD(5); (-t15 * t26 - t17 * t28) * m(6); qJD(4) * t23; m(6) * t27; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

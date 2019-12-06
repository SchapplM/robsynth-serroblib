% Calculate time derivative of joint inertia matrix for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:05
% EndTime: 2019-12-05 15:35:06
% DurationCPUTime: 0.25s
% Computational Cost: add. (166->65), mult. (468->101), div. (0->0), fcn. (325->6), ass. (0->34)
t22 = cos(qJ(4));
t17 = t22 ^ 2;
t20 = sin(qJ(4));
t30 = t20 ^ 2 + t17;
t40 = -mrSges(5,1) - mrSges(6,1);
t39 = mrSges(5,2) - mrSges(6,3);
t12 = -t22 * mrSges(6,1) - t20 * mrSges(6,3);
t38 = -t22 * mrSges(5,1) + t20 * mrSges(5,2) + t12;
t37 = m(6) * qJ(5) + mrSges(6,3);
t36 = -m(6) * pkin(4) + t40;
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t21 = sin(qJ(2));
t23 = cos(qJ(2));
t9 = t18 * t23 + t19 * t21;
t5 = t9 * qJD(2);
t8 = t18 * t21 - t19 * t23;
t35 = t8 * t5;
t6 = t8 * qJD(2);
t34 = t9 * t6;
t13 = t18 * pkin(2) + pkin(6);
t33 = t30 * t13 * t6;
t31 = Ifges(5,4) - Ifges(6,5);
t29 = qJD(4) * t9;
t28 = qJD(4) * t20;
t27 = qJ(5) * qJD(4);
t14 = -t19 * pkin(2) - pkin(3);
t26 = (m(6) * t13 + mrSges(6,2)) * t22;
t24 = -t22 * pkin(4) - t20 * qJ(5);
t11 = (t20 * mrSges(5,1) + t22 * mrSges(5,2)) * qJD(4);
t10 = (t20 * mrSges(6,1) - t22 * mrSges(6,3)) * qJD(4);
t7 = t14 + t24;
t4 = -pkin(4) * t28 + t20 * qJD(5) + t22 * t27;
t1 = [0.2e1 * m(4) * (-t34 + t35) + 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t30 * t34 + t35); (t11 + t10) * t8 + (-t21 * mrSges(3,1) - t23 * mrSges(3,2)) * qJD(2) + (-mrSges(4,1) + t38) * t5 - (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t30) * t6 + m(5) * (t14 * t5 - t33) + m(6) * (-t4 * t8 + t7 * t5 - t33) + m(4) * (-t18 * t6 - t19 * t5) * pkin(2); 0.2e1 * t7 * t10 + 0.2e1 * t14 * t11 + 0.2e1 * (-m(6) * t7 - t12) * t4 + 0.2e1 * t31 * qJD(4) * t17 + 0.2e1 * (-t31 * t20 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t22) * t28; 0; 0; 0; (qJD(5) * t22 - t20 * t27) * t9 * m(6) - (t36 * t20 + (-mrSges(5,2) + t37) * t22) * t6 + (t39 * t20 + t36 * t22) * t29; qJD(5) * t26 + ((-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t22 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t20 + (m(6) * t24 + t38) * t13) * qJD(4); m(6) * t4 + (t40 * t20 - t39 * t22) * qJD(4); 0.2e1 * t37 * qJD(5); (-t20 * t6 + t22 * t29) * m(6); qJD(4) * t26; m(6) * t28; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

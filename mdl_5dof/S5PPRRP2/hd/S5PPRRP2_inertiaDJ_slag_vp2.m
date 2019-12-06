% Calculate time derivative of joint inertia matrix for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:33
% EndTime: 2019-12-05 15:08:34
% DurationCPUTime: 0.22s
% Computational Cost: add. (143->59), mult. (431->91), div. (0->0), fcn. (296->6), ass. (0->33)
t20 = cos(qJ(4));
t15 = t20 ^ 2;
t18 = sin(qJ(4));
t28 = t18 ^ 2 + t15;
t39 = -mrSges(5,1) - mrSges(6,1);
t38 = mrSges(5,2) - mrSges(6,3);
t12 = -t20 * mrSges(6,1) - t18 * mrSges(6,3);
t37 = -t20 * mrSges(5,1) + t18 * mrSges(5,2) + t12;
t36 = m(6) * qJ(5) + mrSges(6,3);
t35 = -m(6) * pkin(4) + t39;
t25 = qJ(5) * qJD(4);
t26 = qJD(4) * t18;
t4 = -pkin(4) * t26 + t18 * qJD(5) + t20 * t25;
t34 = m(6) * t4;
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t8 = t21 * t16 + t19 * t17;
t6 = t8 * qJD(3);
t7 = t19 * t16 - t21 * t17;
t32 = t7 * t6;
t5 = t7 * qJD(3);
t31 = t8 * t5;
t30 = t28 * pkin(6) * t5;
t29 = Ifges(5,4) - Ifges(6,5);
t27 = qJD(4) * t8;
t24 = (m(6) * pkin(6) + mrSges(6,2)) * t20;
t22 = -t20 * pkin(4) - t18 * qJ(5);
t11 = -pkin(3) + t22;
t10 = (t18 * mrSges(5,1) + t20 * mrSges(5,2)) * qJD(4);
t9 = (t18 * mrSges(6,1) - t20 * mrSges(6,3)) * qJD(4);
t1 = [0.2e1 * m(4) * (-t31 + t32) + 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t28 * t31 + t32); 0; 0; (t10 + t9) * t7 + (-mrSges(4,1) + t37) * t6 + m(5) * (-pkin(3) * t6 - t30) + m(6) * (t11 * t6 - t4 * t7 - t30) - (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t28) * t5; 0; -0.2e1 * pkin(3) * t10 - 0.2e1 * t4 * t12 + 0.2e1 * (t9 - t34) * t11 + 0.2e1 * t29 * qJD(4) * t15 + 0.2e1 * (-t29 * t18 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t20) * t26; (qJD(5) * t20 - t18 * t25) * t8 * m(6) - (t35 * t18 + (-mrSges(5,2) + t36) * t20) * t5 + (t38 * t18 + t35 * t20) * t27; t34 + (t39 * t18 - t38 * t20) * qJD(4); qJD(5) * t24 + ((-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t20 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t18 + (m(6) * t22 + t37) * pkin(6)) * qJD(4); 0.2e1 * t36 * qJD(5); (-t18 * t5 + t20 * t27) * m(6); m(6) * t26; qJD(4) * t24; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

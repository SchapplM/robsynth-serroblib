% Calculate time derivative of joint inertia matrix for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:58
% EndTime: 2019-03-08 18:35:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (143->39), mult. (416->64), div. (0->0), fcn. (229->4), ass. (0->30)
t14 = sin(qJ(3));
t32 = pkin(2) * t14;
t16 = cos(qJ(3));
t29 = t16 * pkin(2);
t26 = -mrSges(4,1) - mrSges(5,1);
t31 = m(5) * pkin(3);
t17 = cos(qJ(2));
t13 = t17 * pkin(1) + pkin(2);
t15 = sin(qJ(2));
t27 = t15 * t16;
t10 = pkin(1) * t27 + t14 * t13;
t22 = qJD(3) * t16;
t23 = qJD(3) * t14;
t28 = t14 * t15;
t5 = t13 * t22 + (-t15 * t23 + (t16 * t17 - t28) * qJD(2)) * pkin(1);
t30 = pkin(2) * t10 * t22 + t5 * t32;
t25 = -mrSges(4,2) - mrSges(5,2);
t24 = pkin(2) * qJD(3);
t21 = t25 * t16;
t20 = 0.2e1 * t25;
t9 = -pkin(1) * t28 + t16 * t13;
t6 = -t13 * t23 + (-t15 * t22 + (-t14 * t17 - t27) * qJD(2)) * pkin(1);
t3 = t6 * mrSges(5,1);
t4 = t6 * mrSges(4,1);
t19 = t25 * t5 + t3 + t4;
t18 = (-mrSges(3,1) * t15 - mrSges(3,2) * t17) * qJD(2) * pkin(1);
t12 = pkin(3) + t29;
t8 = pkin(3) + t9;
t1 = t10 * t5;
t2 = [0.2e1 * t3 + 0.2e1 * t4 + t5 * t20 + 0.2e1 * t18 + 0.2e1 * m(4) * (t9 * t6 + t1) + 0.2e1 * m(5) * (t8 * t6 + t1); t18 + (t26 * t14 + t21) * t24 + m(4) * ((t16 * t6 - t9 * t23) * pkin(2) + t30) + m(5) * (-pkin(2) * t8 * t23 + t12 * t6 + t30) + t19; (t20 * t29 + 0.2e1 * ((-t12 + t29) * m(5) + t26) * t32) * qJD(3); t6 * t31 + t19; (t21 + (t26 - t31) * t14) * t24; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1) t2(2) t2(4) t2(7); t2(2) t2(3) t2(5) t2(8); t2(4) t2(5) t2(6) t2(9); t2(7) t2(8) t2(9) t2(10);];
Mq  = res;

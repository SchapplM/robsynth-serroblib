% Calculate time derivative of joint inertia matrix for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:55
% EndTime: 2019-03-08 18:34:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (191->37), mult. (504->64), div. (0->0), fcn. (342->6), ass. (0->31)
t18 = sin(pkin(7));
t23 = cos(qJ(2));
t27 = pkin(1) * qJD(2);
t19 = cos(pkin(7));
t21 = sin(qJ(2));
t28 = t19 * t21;
t10 = (-t18 * t23 - t28) * t27;
t29 = t18 * t21;
t11 = (t19 * t23 - t29) * t27;
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t17 = t23 * pkin(1) + pkin(2);
t13 = pkin(1) * t28 + t18 * t17;
t25 = -pkin(1) * t29 + t19 * t17;
t9 = pkin(3) + t25;
t5 = t22 * t13 + t20 * t9;
t3 = -t5 * qJD(4) + t22 * t10 - t20 * t11;
t1 = t3 * mrSges(5,1);
t34 = t1 + (-mrSges(3,1) * t21 - mrSges(3,2) * t23) * t27 + t10 * mrSges(4,1) - t11 * mrSges(4,2);
t33 = pkin(2) * t18;
t4 = -t20 * t13 + t22 * t9;
t2 = t4 * qJD(4) + t20 * t10 + t22 * t11;
t32 = t2 * mrSges(5,2);
t16 = t19 * pkin(2) + pkin(3);
t14 = t20 * t16 + t22 * t33;
t8 = t14 * qJD(4);
t6 = t8 * mrSges(5,1);
t12 = t22 * t16 - t20 * t33;
t7 = t12 * qJD(4);
t26 = -t7 * mrSges(5,2) - t6;
t15 = [-0.2e1 * t32 + 0.2e1 * m(5) * (t5 * t2 + t4 * t3) + 0.2e1 * m(4) * (t25 * t10 + t13 * t11) + 0.2e1 * t34; -t6 + (-t7 - t2) * mrSges(5,2) + m(5) * (t12 * t3 + t14 * t2 - t8 * t4 + t7 * t5) + m(4) * (t10 * t19 + t11 * t18) * pkin(2) + t34; 0.2e1 * m(5) * (-t12 * t8 + t14 * t7) + 0.2e1 * t26; 0; 0; 0; t1 - t32; t26; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1) t15(2) t15(4) t15(7); t15(2) t15(3) t15(5) t15(8); t15(4) t15(5) t15(6) t15(9); t15(7) t15(8) t15(9) t15(10);];
Mq  = res;

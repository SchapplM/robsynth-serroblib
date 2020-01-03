% Calculate time derivative of joint inertia matrix for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:06
% DurationCPUTime: 0.24s
% Computational Cost: add. (335->60), mult. (749->109), div. (0->0), fcn. (581->6), ass. (0->31)
t36 = qJD(3) + qJD(4);
t35 = 2 * m(5);
t25 = cos(qJ(3));
t28 = -cos(pkin(7)) * pkin(1) - pkin(2);
t34 = -0.2e1 * t25 * pkin(3) + 0.2e1 * t28;
t33 = m(5) * pkin(3);
t20 = sin(pkin(7)) * pkin(1) + pkin(5);
t32 = pkin(6) + t20;
t22 = sin(qJ(4));
t23 = sin(qJ(3));
t24 = cos(qJ(4));
t18 = t22 * t25 + t24 * t23;
t11 = t36 * t18;
t17 = -t22 * t23 + t24 * t25;
t31 = t17 * t11;
t10 = t36 * t17;
t30 = t18 * t10;
t29 = 0.2e1 * t25;
t4 = t11 * mrSges(5,1) + t10 * mrSges(5,2);
t27 = qJD(3) * t32;
t12 = t23 * t27;
t13 = t25 * t27;
t14 = t32 * t23;
t15 = t32 * t25;
t5 = -t24 * t14 - t22 * t15;
t2 = qJD(4) * t5 - t24 * t12 - t22 * t13;
t6 = -t22 * t14 + t24 * t15;
t3 = -qJD(4) * t6 + t22 * t12 - t24 * t13;
t26 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t10 - Ifges(5,6) * t11;
t16 = (-mrSges(5,1) * t22 - mrSges(5,2) * t24) * qJD(4) * pkin(3);
t1 = [0.2e1 * Ifges(5,1) * t30 - 0.2e1 * Ifges(5,2) * t31 + (t6 * t2 + t5 * t3) * t35 + t4 * t34 + ((t28 * mrSges(4,2) + Ifges(4,4) * t25) * t29 + (t33 * t34 + 0.2e1 * pkin(3) * (-t17 * mrSges(5,1) + t18 * mrSges(5,2)) + 0.2e1 * t28 * mrSges(4,1) - 0.2e1 * Ifges(4,4) * t23 + (Ifges(4,1) - Ifges(4,2)) * t29) * t23) * qJD(3) + 0.2e1 * (t10 * t17 - t18 * t11) * Ifges(5,4) + 0.2e1 * (-t5 * t10 - t6 * t11 + t2 * t17 - t3 * t18) * mrSges(5,3); m(5) * (t6 * t10 - t5 * t11 + t3 * t17 + t2 * t18); (t30 - t31) * t35; (Ifges(4,5) * t25 - Ifges(4,6) * t23 + (-mrSges(4,1) * t25 + mrSges(4,2) * t23) * t20) * qJD(3) + (m(5) * (t2 * t22 + t24 * t3 + (-t22 * t5 + t24 * t6) * qJD(4)) + (-t24 * t10 - t22 * t11 + (t17 * t24 + t18 * t22) * qJD(4)) * mrSges(5,3)) * pkin(3) + t26; (-t23 * mrSges(4,1) - t25 * mrSges(4,2)) * qJD(3) + (t10 * t22 - t11 * t24 + (-t17 * t22 + t18 * t24) * qJD(4)) * t33 - t4; 0.2e1 * t16; t26; -t4; t16; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

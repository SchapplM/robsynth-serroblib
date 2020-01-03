% Calculate time derivative of joint inertia matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:22
% DurationCPUTime: 0.30s
% Computational Cost: add. (363->65), mult. (716->111), div. (0->0), fcn. (603->6), ass. (0->35)
t27 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(5));
t32 = cos(qJ(5));
t17 = t32 * t27 + t31 * t29;
t15 = t17 * qJD(5);
t16 = t31 * t27 - t32 * t29;
t49 = t16 * t15;
t28 = sin(pkin(7));
t11 = t16 * t28;
t14 = t16 * qJD(5);
t5 = -t15 * mrSges(6,1) + t14 * mrSges(6,2);
t30 = cos(pkin(7));
t38 = t30 * qJD(2);
t20 = -qJD(4) + t38;
t48 = (t27 ^ 2 + t29 ^ 2) * t20;
t47 = 2 * m(6);
t33 = -pkin(1) - pkin(2);
t40 = t30 * qJ(2) + t28 * t33;
t18 = -qJ(4) + t40;
t44 = pkin(6) - t18;
t39 = qJD(2) * t28;
t35 = -t28 * qJ(2) + t30 * t33;
t34 = pkin(3) - t35;
t8 = t44 * t27;
t9 = t44 * t29;
t3 = t31 * t9 + t32 * t8;
t4 = t31 * t8 - t32 * t9;
t10 = t17 * t28;
t13 = t29 * pkin(4) + t34;
t7 = qJD(5) * t11;
t6 = qJD(5) * t10;
t2 = -t4 * qJD(5) - t17 * t20;
t1 = t3 * qJD(5) - t16 * t20;
t12 = [0.2e1 * t13 * t5 + (t4 * t1 + t13 * t39 + t3 * t2) * t47 + 0.2e1 * Ifges(6,2) * t49 + 0.2e1 * m(5) * (t18 * t48 + t34 * t39) + 0.2e1 * mrSges(4,2) * t38 - 0.2e1 * mrSges(5,3) * t48 + 0.2e1 * (t1 * t16 + t4 * t15 + t2 * t17) * mrSges(6,3) + 0.2e1 * (m(4) * (-t35 * t28 + t40 * t30) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2) + 0.2e1 * (t29 * mrSges(5,1) - t16 * mrSges(6,1) - t27 * mrSges(5,2) - t17 * mrSges(6,2) + mrSges(4,1)) * t39 + 0.2e1 * (-t3 * mrSges(6,3) - Ifges(6,1) * t17) * t14 + 0.2e1 * (t14 * t16 - t17 * t15) * Ifges(6,4); -t30 * t5 + m(5) * (-t38 + t48) * t28 + (t10 * t14 - t11 * t15 - t6 * t16 + t7 * t17) * mrSges(6,3) + (-t11 * t1 - t10 * t2 - t38 * t28 + t7 * t3 - t6 * t4) * m(6); (-t10 * t7 + t11 * t6) * t47; m(6) * (t17 * t1 - t14 * t4 - t15 * t3 - t16 * t2); m(6) * (t15 * t10 + t14 * t11 - t16 * t7 - t17 * t6); (-t17 * t14 + t49) * t47; (m(5) + m(6)) * t39 + t5; 0; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t14 + Ifges(6,6) * t15; t7 * mrSges(6,1) + t6 * mrSges(6,2); t5; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;

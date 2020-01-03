% Calculate time derivative of joint inertia matrix for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:29
% EndTime: 2019-12-31 17:01:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (118->44), mult. (358->71), div. (0->0), fcn. (205->6), ass. (0->33)
t44 = Ifges(5,1) - Ifges(5,2);
t25 = cos(qJ(4));
t33 = qJD(4) * t25;
t23 = sin(qJ(4));
t34 = qJD(4) * t23;
t8 = -mrSges(5,1) * t34 - mrSges(5,2) * t33;
t43 = -0.2e1 * t8;
t40 = Ifges(5,6) * t23;
t21 = sin(pkin(7));
t24 = sin(qJ(2));
t39 = t21 * t24;
t22 = cos(pkin(7));
t38 = t22 * t24;
t26 = cos(qJ(2));
t15 = t26 * pkin(1) + pkin(2);
t37 = pkin(1) * t38 + t21 * t15;
t36 = t23 ^ 2 + t25 ^ 2;
t35 = pkin(1) * qJD(2);
t3 = (t21 * t26 + t38) * t35;
t9 = -t25 * mrSges(5,1) + t23 * mrSges(5,2);
t32 = (-mrSges(4,1) + t9) * t3;
t4 = (t22 * t26 - t39) * t35;
t31 = t36 * t4;
t30 = t36 * mrSges(5,3);
t29 = -pkin(1) * t39 + t22 * t15;
t28 = (-0.2e1 * Ifges(5,4) * t23 + t44 * t25) * t34 + (0.2e1 * Ifges(5,4) * t25 + t44 * t23) * t33;
t27 = (-mrSges(3,1) * t24 - mrSges(3,2) * t26) * t35;
t18 = Ifges(5,5) * t33;
t14 = -t22 * pkin(2) - pkin(3);
t13 = t21 * pkin(2) + pkin(6);
t2 = pkin(6) + t37;
t1 = -pkin(3) - t29;
t5 = [t1 * t43 - 0.2e1 * t4 * mrSges(4,2) + 0.2e1 * m(5) * (t1 * t3 + t2 * t31) + 0.2e1 * m(4) * (-t29 * t3 + t37 * t4) + t28 + 0.2e1 * t30 * t4 + 0.2e1 * t27 + 0.2e1 * t32; (-t1 - t14) * t8 + t32 + t27 + (-mrSges(4,2) + t30) * t4 + m(5) * (t13 * t31 + t14 * t3) + m(4) * (t21 * t4 - t22 * t3) * pkin(2) + t28; t14 * t43 + t28; 0; 0; 0; t18 + (-mrSges(5,1) * t23 - mrSges(5,2) * t25) * t4 + (t9 * t2 - t40) * qJD(4); t18 + (t9 * t13 - t40) * qJD(4); t8; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1), t5(2), t5(4), t5(7); t5(2), t5(3), t5(5), t5(8); t5(4), t5(5), t5(6), t5(9); t5(7), t5(8), t5(9), t5(10);];
Mq = res;

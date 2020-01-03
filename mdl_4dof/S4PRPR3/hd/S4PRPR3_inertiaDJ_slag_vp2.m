% Calculate time derivative of joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:48
% EndTime: 2019-12-31 16:20:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (119->30), mult. (313->57), div. (0->0), fcn. (249->4), ass. (0->20)
t22 = 2 * m(5);
t14 = sin(pkin(7));
t15 = cos(pkin(7));
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t9 = t17 * t14 + t16 * t15;
t7 = t9 * qJD(4);
t8 = -t16 * t14 + t17 * t15;
t21 = t8 * t7;
t6 = t8 * qJD(4);
t20 = t9 * t6;
t19 = pkin(5) + qJ(3);
t10 = t19 * t14;
t11 = t19 * t15;
t4 = -t17 * t10 - t16 * t11;
t5 = -t16 * t10 + t17 * t11;
t3 = t7 * mrSges(5,1) + t6 * mrSges(5,2);
t2 = -t9 * qJD(3) - t5 * qJD(4);
t1 = t8 * qJD(3) + t4 * qJD(4);
t12 = [(t20 - t21) * t22; m(5) * (t1 * t9 + t2 * t8 - t4 * t7 + t5 * t6); 0.2e1 * Ifges(5,1) * t20 - 0.2e1 * Ifges(5,2) * t21 + 0.2e1 * (-t15 * pkin(3) - pkin(2)) * t3 + (t5 * t1 + t4 * t2) * t22 + 0.2e1 * (t6 * t8 - t9 * t7) * Ifges(5,4) + 0.2e1 * (t1 * t8 - t2 * t9 - t4 * t6 - t5 * t7) * mrSges(5,3) + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) * (t14 ^ 2 + t15 ^ 2); 0; t3; 0; -t3; t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t6 - Ifges(5,6) * t7; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t12(1), t12(2), t12(4), t12(7); t12(2), t12(3), t12(5), t12(8); t12(4), t12(5), t12(6), t12(9); t12(7), t12(8), t12(9), t12(10);];
Mq = res;

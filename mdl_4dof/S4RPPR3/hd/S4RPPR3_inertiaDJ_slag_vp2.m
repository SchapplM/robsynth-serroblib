% Calculate time derivative of joint inertia matrix for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (150->32), mult. (344->59), div. (0->0), fcn. (280->6), ass. (0->21)
t23 = 2 * m(5);
t15 = sin(pkin(7));
t16 = cos(pkin(7));
t17 = sin(qJ(4));
t18 = cos(qJ(4));
t10 = -t17 * t15 + t18 * t16;
t11 = t18 * t15 + t17 * t16;
t7 = t11 * qJD(4);
t22 = t10 * t7;
t6 = t10 * qJD(4);
t21 = t11 * t6;
t12 = sin(pkin(6)) * pkin(1) + qJ(3);
t20 = pkin(5) + t12;
t8 = t20 * t15;
t9 = t20 * t16;
t3 = -t17 * t9 - t18 * t8;
t4 = -t17 * t8 + t18 * t9;
t5 = t7 * mrSges(5,1) + t6 * mrSges(5,2);
t2 = -t11 * qJD(3) - t4 * qJD(4);
t1 = t10 * qJD(3) + t3 * qJD(4);
t13 = [(t4 * t1 + t3 * t2) * t23 + 0.2e1 * Ifges(5,1) * t21 - 0.2e1 * Ifges(5,2) * t22 + 0.2e1 * (-cos(pkin(6)) * pkin(1) - pkin(2) - t16 * pkin(3)) * t5 + 0.2e1 * (t6 * t10 - t11 * t7) * Ifges(5,4) + 0.2e1 * (t1 * t10 - t2 * t11 - t3 * t6 - t4 * t7) * mrSges(5,3) + 0.2e1 * (m(4) * t12 + mrSges(4,3)) * qJD(3) * (t15 ^ 2 + t16 ^ 2); m(5) * (t1 * t11 + t2 * t10 - t3 * t7 + t4 * t6); (t21 - t22) * t23; t5; 0; 0; t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t6 - Ifges(5,6) * t7; -t5; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1), t13(2), t13(4), t13(7); t13(2), t13(3), t13(5), t13(8); t13(4), t13(5), t13(6), t13(9); t13(7), t13(8), t13(9), t13(10);];
Mq = res;

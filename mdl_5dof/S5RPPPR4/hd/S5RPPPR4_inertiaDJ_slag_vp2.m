% Calculate time derivative of joint inertia matrix for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:07
% EndTime: 2019-12-31 17:45:07
% DurationCPUTime: 0.17s
% Computational Cost: add. (265->46), mult. (532->76), div. (0->0), fcn. (432->6), ass. (0->28)
t21 = sin(pkin(8));
t23 = cos(pkin(8));
t24 = cos(qJ(5));
t34 = sin(qJ(5));
t11 = -t24 * t21 - t34 * t23;
t26 = t34 * t21 - t24 * t23;
t9 = t26 * qJD(5);
t36 = t11 * t9;
t10 = t11 * qJD(5);
t6 = t26 * t10;
t38 = -t6 + t36;
t29 = (t21 ^ 2 + t23 ^ 2) * qJD(4);
t37 = 2 * m(6);
t15 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t35 = -pkin(6) + t15;
t31 = 2 * mrSges(6,3);
t16 = sin(pkin(7)) * pkin(1) + qJ(3);
t5 = -t9 * mrSges(6,1) + t10 * mrSges(6,2);
t28 = t38 * t37;
t7 = t35 * t21;
t8 = t35 * t23;
t3 = t24 * t8 - t34 * t7;
t27 = -t24 * t7 - t34 * t8;
t1 = t11 * qJD(4) + t3 * qJD(5);
t2 = t26 * qJD(4) + t27 * qJD(5);
t25 = t11 * t1 - t10 * t3 + t2 * t26 - t27 * t9;
t13 = t21 * pkin(4) + t16;
t4 = [-0.2e1 * Ifges(6,1) * t6 + 0.2e1 * Ifges(6,2) * t36 + 0.2e1 * t13 * t5 + (t13 * qJD(3) - t1 * t27 + t3 * t2) * t37 + 0.2e1 * m(5) * (t16 * qJD(3) - t15 * t29) + t25 * t31 + 0.2e1 * (m(4) * t16 + t21 * mrSges(5,1) - t11 * mrSges(6,1) + t23 * mrSges(5,2) - mrSges(6,2) * t26 + mrSges(4,3)) * qJD(3) + 0.2e1 * mrSges(5,3) * t29 + 0.2e1 * (t10 * t11 - t26 * t9) * Ifges(6,4); m(6) * (-t1 * t26 - t10 * t27 + t2 * t11 + t3 * t9); t28; -m(5) * t29 - m(6) * t25 - t38 * t31; 0; t28; (m(5) + m(6)) * qJD(3) + t5; 0; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t9; -t5; t10 * mrSges(6,1) + t9 * mrSges(6,2); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (235->54), mult. (585->84), div. (0->0), fcn. (463->4), ass. (0->27)
t35 = (-mrSges(4,3) - mrSges(5,2));
t21 = cos(pkin(6));
t27 = -t21 * pkin(2) - pkin(1);
t34 = 0.2e1 * t27;
t33 = 2 * t35;
t32 = m(5) * qJ(4) + mrSges(5,3);
t30 = pkin(5) + qJ(2);
t20 = sin(pkin(6));
t15 = t30 * t20;
t16 = t30 * t21;
t22 = sin(qJ(3));
t23 = cos(qJ(3));
t25 = -t23 * t15 - t22 * t16;
t13 = t22 * t20 - t23 * t21;
t3 = -t13 * qJD(2) + t25 * qJD(3);
t14 = t23 * t20 + t22 * t21;
t7 = -t22 * t15 + t23 * t16;
t4 = t14 * qJD(2) + t7 * qJD(3);
t28 = -t25 * t4 + t7 * t3;
t26 = 2 * Ifges(5,5) - 2 * Ifges(4,4);
t11 = t14 * qJD(3);
t10 = t13 * qJD(3);
t9 = t11 * mrSges(5,1);
t8 = t10 * mrSges(4,2);
t5 = t13 * pkin(3) - t14 * qJ(4) + t27;
t2 = t11 * pkin(3) + t10 * qJ(4) - t14 * qJD(4);
t1 = [-t8 * t34 + 0.2e1 * t2 * (t13 * mrSges(5,1) - t14 * mrSges(5,3)) + 0.2e1 * t5 * t9 + 0.2e1 * m(4) * t28 + 0.2e1 * m(5) * (t5 * t2 + t28) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t20 ^ 2 + t21 ^ 2) + (mrSges(4,1) * t34 + t14 * t26 + t7 * t33 + 0.2e1 * (Ifges(5,3) + Ifges(4,2)) * t13) * t11 - (-0.2e1 * t5 * mrSges(5,3) + t13 * t26 + t25 * t33 + 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t14) * t10 - 0.2e1 * t35 * (-t3 * t13 + t4 * t14); m(5) * t2 + t11 * mrSges(4,1) + t10 * mrSges(5,3) - t8 + t9; 0; m(5) * qJD(4) * t7 + (-Ifges(4,6) + Ifges(5,6)) * t11 - (Ifges(4,5) + Ifges(5,4)) * t10 + (pkin(3) * t10 - qJ(4) * t11 - qJD(4) * t13) * mrSges(5,2) + (-m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1)) * t4 + (-mrSges(4,2) + t32) * t3; 0; 0.2e1 * t32 * qJD(4); m(5) * t4 - t10 * mrSges(5,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

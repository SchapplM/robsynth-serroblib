% Calculate time derivative of joint inertia matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:38
% EndTime: 2019-12-31 16:40:38
% DurationCPUTime: 0.21s
% Computational Cost: add. (138->47), mult. (352->84), div. (0->0), fcn. (263->4), ass. (0->22)
t16 = sin(pkin(6));
t17 = cos(pkin(6));
t28 = qJD(2) * (t16 ^ 2 + t17 ^ 2);
t25 = t16 * qJD(3);
t27 = -pkin(5) + qJ(2);
t26 = qJ(2) * t28;
t23 = t16 * qJ(3) + pkin(1);
t18 = sin(qJ(4));
t19 = cos(qJ(4));
t20 = t16 * t18 + t17 * t19;
t5 = t20 * qJD(4);
t9 = t16 * t19 - t17 * t18;
t6 = qJD(4) * t9;
t21 = t6 * mrSges(5,1) - t5 * mrSges(5,2);
t10 = t27 * t16;
t11 = t27 * t17;
t3 = t19 * t10 - t18 * t11;
t4 = t18 * t10 + t19 * t11;
t7 = (pkin(2) + pkin(3)) * t17 + t23;
t2 = t9 * qJD(2) - qJD(4) * t4;
t1 = t20 * qJD(2) + qJD(4) * t3;
t8 = [-0.2e1 * t5 * t9 * Ifges(5,1) + 0.2e1 * t6 * Ifges(5,2) * t20 + 0.2e1 * t7 * t21 + 0.2e1 * (t17 * mrSges(4,1) + mrSges(5,1) * t20 + t9 * mrSges(5,2) + t16 * mrSges(4,3)) * t25 + 0.2e1 * m(5) * (t4 * t1 + t3 * t2 + t7 * t25) + 0.2e1 * m(4) * (-(-t17 * pkin(2) - t23) * t25 + t26) + 0.2e1 * m(3) * t26 + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * t28 + 0.2e1 * (t20 * t5 - t9 * t6) * Ifges(5,4) + 0.2e1 * (-t1 * t20 - t2 * t9 + t3 * t5 - t4 * t6) * mrSges(5,3); (-m(4) - m(5)) * t25 - t21; 0; m(5) * (t18 * t1 + t19 * t2 + (-t18 * t3 + t19 * t4) * qJD(4)) + m(4) * t16 * qJD(2) + (-t18 * t6 + t19 * t5 + (t18 * t9 - t19 * t20) * qJD(4)) * mrSges(5,3); 0; 0; t2 * mrSges(5,1) - t1 * mrSges(5,2) - Ifges(5,5) * t5 - Ifges(5,6) * t6; 0; (-mrSges(5,1) * t18 - mrSges(5,2) * t19) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;

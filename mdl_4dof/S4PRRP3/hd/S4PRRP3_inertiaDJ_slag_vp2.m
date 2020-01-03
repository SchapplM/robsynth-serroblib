% Calculate time derivative of joint inertia matrix for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:44
% DurationCPUTime: 0.12s
% Computational Cost: add. (72->39), mult. (188->59), div. (0->0), fcn. (96->2), ass. (0->17)
t19 = -2 * pkin(2);
t9 = cos(qJ(3));
t18 = -0.2e1 * t9 * pkin(3) - (2 * pkin(2));
t17 = -2 * mrSges(5,3);
t16 = m(5) * pkin(3);
t8 = sin(qJ(3));
t12 = qJD(3) * t8;
t15 = qJD(3) * t9 * mrSges(5,2) + mrSges(5,1) * t12;
t14 = Ifges(5,4) + Ifges(4,4);
t13 = -qJ(4) - pkin(5);
t11 = 0.2e1 * t9;
t10 = qJD(3) * t13;
t4 = t13 * t9;
t3 = t13 * t8;
t2 = -t8 * qJD(4) + t9 * t10;
t1 = t9 * qJD(4) + t8 * t10;
t5 = [0; m(5) * (t1 * t8 + t2 * t9 + (-t3 * t8 - t4 * t9) * qJD(3)); t15 * t18 + 0.2e1 * m(5) * (-t4 * t1 + t3 * t2) + 0.2e1 * (t1 * t9 - t2 * t8) * mrSges(5,3) + (((mrSges(4,2) * t19) + t14 * t11 + t3 * t17) * t9 + (-t4 * t17 + (mrSges(4,1) * t19) + t16 * t18 + 0.2e1 * (pkin(3) * mrSges(5,2) - t14) * t8 + (-pkin(3) * mrSges(5,1) + Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,2)) * t11) * t8) * qJD(3); (-mrSges(4,2) * t9 + (-mrSges(4,1) - t16) * t8) * qJD(3) - t15; -t1 * mrSges(5,2) + (mrSges(5,1) + t16) * t2 + ((mrSges(4,2) * pkin(5) - Ifges(4,6) - Ifges(5,6)) * t8 + (-(mrSges(4,1) * pkin(5)) - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t9) * qJD(3); 0; 0; t12 * t16 + t15; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1), t5(2), t5(4), t5(7); t5(2), t5(3), t5(5), t5(8); t5(4), t5(5), t5(6), t5(9); t5(7), t5(8), t5(9), t5(10);];
Mq = res;

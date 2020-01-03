% Calculate time derivative of joint inertia matrix for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:45
% DurationCPUTime: 0.19s
% Computational Cost: add. (140->54), mult. (266->78), div. (0->0), fcn. (143->4), ass. (0->23)
t23 = 2 * mrSges(6,3);
t12 = cos(qJ(4));
t15 = 0.2e1 * t12;
t22 = 2 * mrSges(5,1);
t21 = (m(6) * pkin(4));
t20 = mrSges(5,2) + mrSges(6,2);
t19 = Ifges(6,4) + Ifges(5,4);
t7 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t18 = qJ(5) - t7;
t11 = sin(qJ(4));
t17 = qJD(4) * t11;
t16 = qJD(4) * t12;
t14 = mrSges(6,1) + t21;
t4 = t18 * t12;
t8 = sin(pkin(7)) * pkin(1) + qJ(3);
t1 = -t12 * qJD(5) + t18 * t17;
t2 = -qJD(4) * t4 - t11 * qJD(5);
t13 = -t12 * t1 - t11 * t2;
t9 = mrSges(6,1) * t16;
t6 = pkin(4) * t16 + qJD(3);
t5 = t11 * pkin(4) + t8;
t3 = t18 * t11;
t10 = [0.2e1 * t6 * (t11 * mrSges(6,1) + t12 * mrSges(6,2)) + 0.2e1 * t5 * t9 + 0.2e1 * m(6) * (-t4 * t1 - t3 * t2 + t5 * t6) + t13 * t23 + (t11 * t22 + mrSges(5,2) * t15 + (2 * mrSges(4,3)) + 0.2e1 * (m(4) + m(5)) * t8) * qJD(3) + ((-t19 * t15 + t8 * t22 + t3 * t23) * t12 + (-0.2e1 * t8 * mrSges(5,2) - 0.2e1 * t5 * mrSges(6,2) - 0.2e1 * t4 * mrSges(6,3) + 0.2e1 * t19 * t11 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t15) * t11) * qJD(4); m(6) * (-t1 * t11 + t2 * t12 + (t11 * t3 + t12 * t4) * qJD(4)); 0; m(6) * ((t4 * t11 - t3 * t12) * qJD(4) - t13); 0; 0; -t2 * mrSges(6,2) + t14 * t1 + ((-mrSges(5,2) * t7 - Ifges(5,6) - Ifges(6,6)) * t12 + (-mrSges(5,1) * t7 + mrSges(6,3) * pkin(4) - Ifges(5,5) - Ifges(6,5)) * t11) * qJD(4); -t9 + ((-mrSges(5,1) - t21) * t12 + t20 * t11) * qJD(4); (-t20 * t12 + (-mrSges(5,1) - t14) * t11) * qJD(4); 0; m(6) * t6 - mrSges(6,2) * t17 + t9; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;

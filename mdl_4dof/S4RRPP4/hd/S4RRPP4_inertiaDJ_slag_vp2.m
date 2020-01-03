% Calculate time derivative of joint inertia matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:58:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (120->64), mult. (296->83), div. (0->0), fcn. (139->2), ass. (0->26)
t27 = 2 * mrSges(5,3);
t26 = -2 * pkin(1);
t11 = cos(qJ(2));
t10 = sin(qJ(2));
t20 = t10 * qJ(3);
t14 = -t11 * pkin(2) - t20;
t6 = -pkin(1) + t14;
t25 = 0.2e1 * t6;
t24 = pkin(2) + pkin(3);
t23 = -mrSges(4,2) + mrSges(5,3);
t22 = pkin(5) - qJ(4);
t21 = qJ(3) * t11;
t19 = qJD(2) * t10;
t18 = t10 * qJD(3);
t17 = 0.2e1 * t11;
t16 = Ifges(3,4) - Ifges(5,4) - Ifges(4,5);
t8 = t22 * t11;
t15 = -t11 * mrSges(4,1) - t10 * mrSges(4,3);
t13 = (m(4) * pkin(5) - t23) * t11;
t9 = qJD(2) * t11 * mrSges(5,2);
t7 = t22 * t10;
t5 = t24 * t11 + pkin(1) + t20;
t4 = qJD(2) * t8 - t10 * qJD(4);
t2 = -t11 * qJD(4) - t22 * t19;
t1 = t18 + (-t24 * t10 + t21) * qJD(2);
t3 = [0.2e1 * t1 * (t11 * mrSges(5,1) + t10 * mrSges(5,2)) + 0.2e1 * t5 * t9 + 0.2e1 * m(5) * (t5 * t1 + t8 * t2 + t7 * t4) + (-t4 * t10 - t2 * t11) * t27 + (((mrSges(3,2) * t26) - 0.2e1 * t6 * mrSges(4,3) - 0.2e1 * t7 * mrSges(5,3) + t16 * t17) * t11 + ((mrSges(3,1) * t26) + mrSges(4,1) * t25 - 0.2e1 * t5 * mrSges(5,1) + t8 * t27 - 0.2e1 * t16 * t10 + (Ifges(3,1) + Ifges(4,1) + Ifges(5,1) - Ifges(3,2) - Ifges(5,2) - Ifges(4,3)) * t17) * t10) * qJD(2) + (m(4) * t25 + 0.2e1 * t15) * (-t18 + (pkin(2) * t10 - t21) * qJD(2)); m(5) * (qJ(3) * t2 + qJD(3) * t8 - t24 * t4) - t4 * mrSges(5,1) + t2 * mrSges(5,2) + qJD(3) * t13 + ((-pkin(2) * mrSges(4,2) + mrSges(5,3) * t24 + Ifges(4,4) + Ifges(3,5) - Ifges(5,5)) * t11 + (t23 * qJ(3) - Ifges(3,6) + Ifges(4,6) - Ifges(5,6)) * t10 + (m(4) * t14 - t11 * mrSges(3,1) + t10 * mrSges(3,2) + t15) * pkin(5)) * qJD(2); 0.2e1 * (mrSges(5,2) + mrSges(4,3) + (m(4) + m(5)) * qJ(3)) * qJD(3); m(5) * t4 + qJD(2) * t13; 0; 0; m(5) * t1 - mrSges(5,1) * t19 + t9; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;

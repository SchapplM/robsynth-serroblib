% Calculate time derivative of joint inertia matrix for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:47
% EndTime: 2019-12-31 17:36:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (176->55), mult. (462->96), div. (0->0), fcn. (355->4), ass. (0->25)
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t32 = qJD(3) * (t18 ^ 2 + t19 ^ 2);
t26 = t18 * qJD(4);
t31 = 2 * m(6);
t20 = sin(qJ(5));
t21 = cos(qJ(5));
t22 = t18 * t20 + t19 * t21;
t11 = t18 * t21 - t19 * t20;
t8 = qJD(5) * t11;
t30 = t22 * t8;
t7 = t22 * qJD(5);
t29 = t11 * t7;
t28 = -pkin(6) + qJ(3);
t27 = qJ(3) * t32;
t3 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t24 = t18 * qJ(4) + pkin(2);
t12 = t28 * t18;
t13 = t28 * t19;
t4 = t21 * t12 - t20 * t13;
t5 = t20 * t12 + t21 * t13;
t9 = (pkin(3) + pkin(4)) * t19 + t24;
t2 = t11 * qJD(3) - qJD(5) * t5;
t1 = t22 * qJD(3) + qJD(5) * t4;
t6 = [(-t29 + t30) * t31; m(6) * (t1 * t11 - t2 * t22 - t4 * t8 - t5 * t7); -0.2e1 * Ifges(6,1) * t29 + 0.2e1 * Ifges(6,2) * t30 + 0.2e1 * t9 * t3 + 0.2e1 * (t19 * mrSges(5,1) + mrSges(6,1) * t22 + t11 * mrSges(6,2) + t18 * mrSges(5,3)) * t26 + (t5 * t1 + t4 * t2 + t9 * t26) * t31 + 0.2e1 * m(4) * t27 + 0.2e1 * m(5) * (-(-t19 * pkin(3) - t24) * t26 + t27) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * t32 + 0.2e1 * (-t8 * t11 + t22 * t7) * Ifges(6,4) + 0.2e1 * (-t1 * t22 - t2 * t11 + t4 * t7 - t5 * t8) * mrSges(6,3); 0; (-m(5) - m(6)) * t26 - t3; 0; m(6) * (-t20 * t7 - t21 * t8 + (t11 * t21 + t20 * t22) * qJD(5)); m(6) * (t20 * t1 + t21 * t2 + (-t20 * t4 + t21 * t5) * qJD(5)) + m(5) * t18 * qJD(3) + (-t20 * t8 + t21 * t7 + (t11 * t20 - t21 * t22) * qJD(5)) * mrSges(6,3); 0; 0; -t3; t2 * mrSges(6,1) - t1 * mrSges(6,2) - Ifges(6,5) * t7 - Ifges(6,6) * t8; 0; (-mrSges(6,1) * t20 - mrSges(6,2) * t21) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;

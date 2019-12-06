% Calculate time derivative of joint inertia matrix for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:00
% EndTime: 2019-12-05 15:03:01
% DurationCPUTime: 0.13s
% Computational Cost: add. (86->36), mult. (260->61), div. (0->0), fcn. (184->6), ass. (0->22)
t13 = sin(qJ(5));
t15 = cos(qJ(5));
t7 = t13 * mrSges(6,1) + t15 * mrSges(6,2);
t28 = -mrSges(5,3) - t7;
t12 = sin(pkin(8));
t14 = sin(qJ(3));
t22 = cos(pkin(8));
t25 = cos(qJ(3));
t4 = t14 * t12 - t25 * t22;
t11 = t15 ^ 2;
t2 = t4 * qJD(3);
t5 = t25 * t12 + t14 * t22;
t1 = t5 * t2;
t23 = t13 ^ 2 + t11;
t21 = qJD(5) * t13;
t20 = qJD(5) * t15;
t3 = t5 * qJD(3);
t19 = t23 * t3;
t17 = -qJ(4) * t2 + t5 * qJD(4);
t16 = -pkin(3) - pkin(6);
t6 = -mrSges(6,1) * t20 + mrSges(6,2) * t21;
t8 = [0.2e1 * m(6) * (t4 * t19 - t1) + 0.2e1 * (m(4) + m(5)) * (t4 * t3 - t1); 0; 0; -t5 * t6 + (mrSges(4,2) + t28) * t2 + (-t23 * mrSges(6,3) - mrSges(4,1) + mrSges(5,2)) * t3 + m(5) * (-pkin(3) * t3 + t17) + m(6) * (t16 * t19 + t17); 0; -0.2e1 * t11 * Ifges(6,4) * qJD(5) - 0.2e1 * qJ(4) * t6 + 0.2e1 * ((m(5) + m(6)) * qJ(4) - t28) * qJD(4) + 0.2e1 * (Ifges(6,4) * t13 + (-Ifges(6,1) + Ifges(6,2)) * t15) * t21; 0.2e1 * (m(5) / 0.2e1 + m(6) * t23 / 0.2e1) * t3; 0; 0; 0; (-t13 * t3 - t4 * t20) * mrSges(6,2) + (t15 * t3 - t4 * t21) * mrSges(6,1); t6; ((-mrSges(6,2) * t16 - Ifges(6,6)) * t15 + (-mrSges(6,1) * t16 - Ifges(6,5)) * t13) * qJD(5); -t7 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;

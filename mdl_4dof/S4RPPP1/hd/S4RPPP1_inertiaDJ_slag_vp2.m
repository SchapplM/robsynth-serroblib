% Calculate time derivative of joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:19
% EndTime: 2019-03-08 18:26:20
% DurationCPUTime: 0.22s
% Computational Cost: add. (54->41), mult. (179->70), div. (0->0), fcn. (127->4), ass. (0->23)
t24 = m(4) + m(5);
t7 = sin(pkin(6));
t17 = t7 * qJD(3);
t8 = sin(pkin(4));
t14 = t8 * t17;
t18 = qJ(2) * t8;
t10 = cos(pkin(4));
t19 = pkin(1) * t10;
t9 = cos(pkin(6));
t23 = t9 * t18 + t7 * t19;
t22 = mrSges(4,1) * t8;
t21 = mrSges(5,1) * t8;
t20 = mrSges(3,3) * t8;
t16 = t8 * qJD(2);
t15 = t7 * t16;
t13 = t9 * t16;
t12 = -pkin(1) * t9 - pkin(2);
t11 = -qJ(3) * t7 - pkin(1);
t4 = t7 * t18;
t3 = t10 * qJD(3) + t13;
t2 = -t10 * qJD(4) + t15;
t1 = (-qJD(4) * t9 - t17) * t8;
t5 = [0.2e1 * m(3) * (t23 * t9 - (t9 * t19 - t4) * t7) * t16 + 0.2e1 * (-t10 * mrSges(3,2) + t9 * t20) * t13 + 0.2e1 * (m(5) * (t4 + (-qJ(4) + t12) * t10) - t10 * mrSges(5,3) + t7 * t21) * t2 + 0.2e1 * ((t22 + t21) * t9 + t24 * t23 + (t24 * qJ(3) + mrSges(5,2) + mrSges(4,3)) * t10) * t3 + 0.2e1 * (m(4) * ((t12 * t10 + t4) * t7 * qJD(2) - (-pkin(2) * t9 + t11) * t14) + t1 * (-mrSges(5,2) * t7 - mrSges(5,3) * t9) + m(5) * (((-pkin(2) - qJ(4)) * t9 + t11) * t1 + (t7 * t2 + t9 * t3) * pkin(3)) - (mrSges(4,2) * t9 - mrSges(4,3) * t7) * t14) * t8 + 0.2e1 * ((t20 + t22) * t7 + (-mrSges(3,1) + mrSges(4,2)) * t10) * t15; -m(4) * t14 + m(5) * t1; 0; m(4) * t15 + m(5) * t2; 0; 0; m(5) * t3; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;

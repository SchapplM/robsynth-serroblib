% Calculate time derivative of joint inertia matrix for
% S4RRPR2
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:31
% EndTime: 2019-07-18 18:16:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (206->54), mult. (376->79), div. (0->0), fcn. (198->4), ass. (0->27)
t15 = sin(qJ(2));
t17 = cos(qJ(2));
t22 = pkin(1) * qJD(2);
t30 = (-mrSges(3,2) * t17 + (-mrSges(3,1) - mrSges(4,1)) * t15) * t22;
t29 = 2 * m(5);
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t21 = t15 * t22;
t20 = -t17 * pkin(1) - pkin(2);
t10 = -pkin(3) + t20;
t11 = t15 * pkin(1) + qJ(3);
t3 = t16 * t10 - t14 * t11;
t9 = t17 * t22 + qJD(3);
t1 = qJD(4) * t3 + t14 * t21 + t16 * t9;
t28 = t1 * mrSges(5,2);
t4 = t14 * t10 + t16 * t11;
t2 = -qJD(4) * t4 - t14 * t9 + t16 * t21;
t27 = t2 * mrSges(5,1);
t18 = -pkin(2) - pkin(3);
t7 = -t14 * qJ(3) + t16 * t18;
t5 = t16 * qJD(3) + qJD(4) * t7;
t26 = t5 * mrSges(5,2);
t8 = t16 * qJ(3) + t14 * t18;
t6 = -t14 * qJD(3) - qJD(4) * t8;
t25 = t6 * mrSges(5,1);
t23 = (-mrSges(5,1) * t14 - mrSges(5,2) * t16) * qJD(4);
t12 = [-0.2e1 * t27 + 0.2e1 * t28 + 0.2e1 * t9 * mrSges(4,3) + 0.2e1 * t30 + (t4 * t1 + t3 * t2) * t29 + 0.2e1 * m(4) * (t11 * t9 + t20 * t21); (qJD(3) + t9) * mrSges(4,3) + (t5 + t1) * mrSges(5,2) + (-t2 - t6) * mrSges(5,1) + t30 + m(5) * (t8 * t1 + t7 * t2 + t6 * t3 + t5 * t4) + m(4) * (-pkin(2) * t21 + qJ(3) * t9 + qJD(3) * t11); 0.2e1 * t26 - 0.2e1 * t25 + (t8 * t5 + t7 * t6) * t29 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3); m(5) * (t14 * t1 + t16 * t2 + (-t14 * t3 + t16 * t4) * qJD(4)) + m(4) * t21 - t23; m(5) * (t14 * t5 + t16 * t6 + (-t14 * t7 + t16 * t8) * qJD(4)) - t23; 0; t27 - t28; t25 - t26; t23; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t12(1), t12(2), t12(4), t12(7); t12(2), t12(3), t12(5), t12(8); t12(4), t12(5), t12(6), t12(9); t12(7), t12(8), t12(9), t12(10);];
Mq  = res;

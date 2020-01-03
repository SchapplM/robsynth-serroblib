% Calculate time derivative of joint inertia matrix for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:39
% EndTime: 2019-12-31 17:35:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (191->51), mult. (549->92), div. (0->0), fcn. (385->6), ass. (0->33)
t22 = sin(qJ(5));
t25 = cos(qJ(5));
t36 = t22 ^ 2 + t25 ^ 2;
t23 = sin(qJ(4));
t24 = sin(qJ(3));
t26 = cos(qJ(4));
t27 = cos(qJ(3));
t11 = t23 * t24 - t26 * t27;
t51 = qJD(3) + qJD(4);
t5 = t51 * t11;
t32 = t36 * t5;
t14 = -t25 * mrSges(6,1) + t22 * mrSges(6,2);
t55 = -mrSges(5,1) + t14;
t54 = Ifges(6,1) - Ifges(6,2);
t52 = t36 * t26;
t50 = (mrSges(6,3) * t36 - mrSges(5,2)) * t26 + t55 * t23;
t49 = m(5) / 0.2e1;
t12 = t23 * t27 + t26 * t24;
t6 = t51 * t12;
t48 = t11 * t6;
t43 = Ifges(6,6) * t22;
t42 = t11 * t23;
t35 = pkin(3) * qJD(4);
t34 = qJD(5) * t22;
t33 = qJD(5) * t25;
t30 = mrSges(6,1) * t22 + mrSges(6,2) * t25;
t29 = (-0.2e1 * Ifges(6,4) * t22 + t25 * t54) * t34 + (0.2e1 * Ifges(6,4) * t25 + t22 * t54) * t33;
t13 = t30 * qJD(5);
t28 = t5 * mrSges(5,2) - mrSges(6,3) * t32 + t11 * t13 + t55 * t6;
t19 = Ifges(6,5) * t33;
t18 = -t26 * pkin(3) - pkin(4);
t17 = t23 * pkin(3) + pkin(7);
t1 = [0; 0; 0.2e1 * m(5) * (-t12 * t5 + t48) + 0.2e1 * m(6) * (-t12 * t32 + t48); 0; m(6) * (-t17 * t32 + t18 * t6) + (-t24 * mrSges(4,1) - t27 * mrSges(4,2)) * qJD(3) + 0.2e1 * ((-t23 * t5 - t26 * t6) * t49 + ((t12 * t26 + t42) * t49 + m(6) * (t52 * t12 + t42) / 0.2e1) * qJD(4)) * pkin(3) + t28; 0.2e1 * t18 * t13 + 0.2e1 * (m(6) * (t52 * t17 + t18 * t23) + t50) * t35 + t29; 0; m(6) * (-pkin(4) * t6 - pkin(7) * t32) + t28; (-pkin(4) + t18) * t13 + (m(6) * (-pkin(4) * t23 + t52 * pkin(7)) + t50) * t35 + t29; -0.2e1 * pkin(4) * t13 + t29; t13; (t12 * t34 + t25 * t5) * mrSges(6,2) + (-t12 * t33 + t22 * t5) * mrSges(6,1); t19 - t30 * t26 * t35 + (t14 * t17 - t43) * qJD(5); t19 + (t14 * pkin(7) - t43) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:28
% DurationCPUTime: 0.36s
% Computational Cost: add. (221->68), mult. (784->134), div. (0->0), fcn. (693->8), ass. (0->43)
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t38 = cos(pkin(8));
t19 = sin(pkin(8));
t20 = cos(pkin(9));
t41 = t19 * t20;
t10 = -t38 * t22 + t24 * t41;
t8 = qJD(4) * t10;
t51 = m(6) * pkin(6);
t23 = cos(qJ(5));
t17 = t23 ^ 2;
t21 = sin(qJ(5));
t40 = t21 ^ 2 + t17;
t50 = mrSges(6,3) + t51;
t14 = -t23 * mrSges(6,1) + t21 * mrSges(6,2);
t49 = -m(6) * pkin(4) - mrSges(5,1) + t14;
t48 = 0.2e1 * m(6);
t46 = m(6) / 0.2e1;
t9 = t22 * t41 + t38 * t24;
t45 = t9 * t8;
t18 = sin(pkin(9));
t43 = t18 * t19;
t42 = t18 * t24;
t39 = qJD(4) * t9;
t37 = qJD(4) * t22;
t36 = qJD(4) * t24;
t35 = qJD(5) * t21;
t34 = qJD(5) * t23;
t32 = t18 * t37;
t30 = t22 * t36;
t4 = t10 * t23 + t21 * t43;
t1 = -qJD(5) * t4 + t21 * t39;
t3 = -t10 * t21 + t23 * t43;
t2 = qJD(5) * t3 - t23 * t39;
t28 = -t1 * t21 + t2 * t23;
t27 = t22 * t8 + t9 * t36;
t12 = -t21 * t20 + t23 * t42;
t11 = -t23 * t20 - t21 * t42;
t5 = qJD(5) * t11 - t23 * t32;
t6 = -qJD(5) * t12 + t21 * t32;
t25 = -t6 * t21 + t5 * t23 + (-t11 * t23 - t12 * t21) * qJD(5);
t13 = (mrSges(6,1) * t21 + mrSges(6,2) * t23) * qJD(5);
t7 = [0.2e1 * m(5) * (-t10 * t39 + t45) + 0.2e1 * m(6) * (t3 * t1 + t4 * t2 + t45); m(6) * (t11 * t1 + t12 * t2 + t6 * t3 + t5 * t4) + 0.2e1 * ((-t10 * t37 - t24 * t39 + t27) * m(5) / 0.2e1 + t27 * t46) * t18; (t18 ^ 2 * t30 + t11 * t6 + t12 * t5) * t48; 0.2e1 * ((-t8 + (-t21 * t3 + t23 * t4) * qJD(4)) * t24 + (-t3 * t34 - t4 * t35 + t28 + t39) * t22) * t46; m(6) * ((-t11 * t21 + t12 * t23 - t42) * t36 + (t25 + t32) * t22); (-0.1e1 + t40) * t30 * t48; t39 * mrSges(5,2) + t9 * t13 + t50 * ((-t21 * t4 - t23 * t3) * qJD(5) + t28) + t49 * t8; t50 * t25 + ((qJD(4) * mrSges(5,2) + t13) * t22 + t49 * t36) * t18; (t40 * t24 * t51 + t49 * t22) * qJD(4) + (-t13 + (t40 * mrSges(6,3) - mrSges(5,2)) * qJD(4)) * t24; 0.2e1 * qJD(5) * t17 * Ifges(6,4) - 0.2e1 * pkin(4) * t13 + 0.2e1 * (-Ifges(6,4) * t21 + (Ifges(6,1) - Ifges(6,2)) * t23) * t35; t1 * mrSges(6,1) - t2 * mrSges(6,2); t6 * mrSges(6,1) - t5 * mrSges(6,2); (t22 * t35 - t23 * t36) * mrSges(6,2) + (-t21 * t36 - t22 * t34) * mrSges(6,1); (Ifges(6,5) * t23 - Ifges(6,6) * t21 + t14 * pkin(6)) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;

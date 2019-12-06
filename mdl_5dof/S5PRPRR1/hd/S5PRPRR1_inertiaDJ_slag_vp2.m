% Calculate time derivative of joint inertia matrix for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:40
% EndTime: 2019-12-05 15:42:41
% DurationCPUTime: 0.40s
% Computational Cost: add. (842->91), mult. (1935->153), div. (0->0), fcn. (1827->6), ass. (0->46)
t38 = sin(pkin(9));
t48 = pkin(6) + qJ(3);
t33 = t48 * t38;
t39 = cos(pkin(9));
t34 = t48 * t39;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t21 = -t41 * t33 + t43 * t34;
t56 = m(6) * pkin(4);
t32 = t43 * t38 + t41 * t39;
t25 = t32 * qJD(4);
t55 = t25 * pkin(4);
t49 = t43 * t39;
t31 = -t41 * t38 + t49;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t19 = t40 * t31 + t42 * t32;
t24 = t31 * qJD(4);
t11 = -qJD(5) * t19 - t40 * t24 - t42 * t25;
t18 = t42 * t31 - t40 * t32;
t54 = t18 * t11;
t10 = qJD(5) * t18 + t42 * t24 - t40 * t25;
t53 = t19 * t10;
t52 = t25 * mrSges(5,1);
t51 = t31 * t25;
t50 = t32 * t24;
t29 = t43 * t33;
t46 = -t39 * pkin(3) - pkin(2);
t4 = -t11 * mrSges(6,1) + t10 * mrSges(6,2);
t20 = -t41 * t34 - t29;
t23 = t24 * mrSges(5,2);
t45 = t23 + t4;
t14 = -qJD(4) * t29 + qJD(3) * t49 + (-qJD(3) * t38 - qJD(4) * t34) * t41;
t12 = -t25 * pkin(7) + t14;
t15 = -t32 * qJD(3) - t21 * qJD(4);
t13 = -t24 * pkin(7) + t15;
t16 = -t32 * pkin(7) + t20;
t17 = t31 * pkin(7) + t21;
t5 = t42 * t16 - t40 * t17;
t2 = qJD(5) * t5 + t42 * t12 + t40 * t13;
t6 = t40 * t16 + t42 * t17;
t3 = -qJD(5) * t6 - t40 * t12 + t42 * t13;
t44 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t11;
t27 = (-mrSges(6,1) * t40 - mrSges(6,2) * t42) * qJD(5) * pkin(4);
t22 = -t31 * pkin(4) + t46;
t1 = [0.2e1 * m(5) * (t50 - t51) + 0.2e1 * m(6) * (t53 + t54); m(5) * (t14 * t32 + t15 * t31 - t20 * t25 + t21 * t24) + m(6) * (t6 * t10 + t5 * t11 + t3 * t18 + t2 * t19); 0.2e1 * Ifges(6,1) * t53 + 0.2e1 * Ifges(6,2) * t54 + 0.2e1 * Ifges(5,1) * t50 - 0.2e1 * Ifges(5,2) * t51 + 0.2e1 * (-t18 * mrSges(6,1) + t19 * mrSges(6,2)) * t55 + 0.2e1 * t22 * t4 + 0.2e1 * t46 * (t23 + t52) + 0.2e1 * m(6) * (t6 * t2 + t22 * t55 + t5 * t3) + 0.2e1 * m(5) * (t21 * t14 + t20 * t15) + 0.2e1 * (t10 * t18 + t19 * t11) * Ifges(6,4) + 0.2e1 * (t24 * t31 - t32 * t25) * Ifges(5,4) + 0.2e1 * (-t5 * t10 + t6 * t11 + t2 * t18 - t3 * t19) * mrSges(6,3) + 0.2e1 * (t14 * t31 - t15 * t32 - t20 * t24 - t21 * t25) * mrSges(5,3) + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) * (t38 ^ 2 + t39 ^ 2); 0; -(-mrSges(5,1) - t56) * t25 + t45; 0; -t52 + (t10 * t40 + t11 * t42 + (-t18 * t40 + t19 * t42) * qJD(5)) * t56 - t45; t15 * mrSges(5,1) - t14 * mrSges(5,2) + Ifges(5,5) * t24 - Ifges(5,6) * t25 + (m(6) * (t2 * t40 + t3 * t42 + (-t40 * t5 + t42 * t6) * qJD(5)) + (-t42 * t10 + t40 * t11 + (t18 * t42 + t19 * t40) * qJD(5)) * mrSges(6,3)) * pkin(4) + t44; 0; 0.2e1 * t27; -t4; t44; 0; t27; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

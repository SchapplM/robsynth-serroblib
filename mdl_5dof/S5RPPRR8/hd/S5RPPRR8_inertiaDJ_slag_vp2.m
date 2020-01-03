% Calculate time derivative of joint inertia matrix for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:05
% EndTime: 2019-12-31 18:01:05
% DurationCPUTime: 0.30s
% Computational Cost: add. (443->76), mult. (826->116), div. (0->0), fcn. (602->6), ass. (0->44)
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t11 = t31 * t28 - t33 * t29;
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t37 = mrSges(6,1) * t30 + mrSges(6,2) * t32;
t14 = t37 * qJD(5);
t46 = t30 ^ 2 + t32 ^ 2;
t40 = mrSges(6,3) * t46;
t7 = t11 * qJD(4);
t59 = t11 * t14 + (mrSges(5,2) - t40) * t7;
t34 = -pkin(1) - pkin(2);
t38 = -t28 * qJ(2) + t29 * t34;
t13 = -pkin(3) + t38;
t17 = t29 * qJ(2) + t28 * t34;
t5 = t33 * t13 - t31 * t17;
t1 = -t11 * qJD(2) + t5 * qJD(4);
t57 = t1 * t40;
t18 = -t32 * mrSges(6,1) + t30 * mrSges(6,2);
t56 = mrSges(5,1) - t18;
t19 = Ifges(6,4) * t30 + Ifges(6,2) * t32;
t20 = Ifges(6,1) * t30 + Ifges(6,4) * t32;
t44 = qJD(5) * t32;
t45 = qJD(5) * t30;
t36 = t32 * (Ifges(6,4) * t44 - Ifges(6,2) * t45) + t30 * (Ifges(6,1) * t44 - Ifges(6,4) * t45);
t55 = (t30 * t19 - t32 * t20) * qJD(5) - t36;
t54 = -0.2e1 * t14;
t52 = t1 * mrSges(5,2);
t12 = t33 * t28 + t31 * t29;
t6 = t31 * t13 + t33 * t17;
t2 = t12 * qJD(2) + t6 * qJD(4);
t51 = t11 * t2;
t8 = t12 * qJD(4);
t50 = t11 * t8;
t49 = t2 * mrSges(5,1);
t48 = t2 * t18;
t42 = t46 * t1;
t41 = t46 * t7;
t39 = t46 * t12;
t4 = -pkin(7) + t6;
t3 = pkin(4) - t5;
t9 = [0.2e1 * m(6) * (t3 * t2 + t4 * t42) - 0.2e1 * t48 + t3 * t54 + t20 * t44 - t19 * t45 + 0.2e1 * m(5) * (t6 * t1 - t5 * t2) + 0.2e1 * t52 + 0.2e1 * t49 + t36 - 0.2e1 * t57 + 0.2e1 * (m(4) * (t17 * t29 - t38 * t28) + m(3) * qJ(2) + mrSges(3,3) + t29 * mrSges(4,2) + t28 * mrSges(4,1)) * qJD(2); t56 * t8 + m(6) * (t1 * t39 + t8 * t3 - t4 * t41 + t51) + m(5) * (t12 * t1 - t8 * t5 - t7 * t6 + t51) - t59; 0.2e1 * m(5) * (-t12 * t7 + t50) + 0.2e1 * m(6) * (-t7 * t39 + t50); 0; 0; 0; m(6) * (-pkin(4) * t2 + pkin(7) * t42) + t48 - t52 - t49 + (t3 + pkin(4)) * t14 + t57 + t55; -m(6) * pkin(7) * t41 + (-m(6) * pkin(4) - t56) * t8 + t59; 0; pkin(4) * t54 - t55; -t37 * t1 + ((-mrSges(6,1) * t4 - Ifges(6,5)) * t32 + (mrSges(6,2) * t4 + Ifges(6,6)) * t30) * qJD(5); (t12 * t45 + t32 * t7) * mrSges(6,2) + (-t12 * t44 + t30 * t7) * mrSges(6,1); -t14; (Ifges(6,5) * t32 - Ifges(6,6) * t30 + t18 * pkin(7)) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;

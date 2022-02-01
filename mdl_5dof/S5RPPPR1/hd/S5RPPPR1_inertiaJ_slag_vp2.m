% Calculate joint inertia matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:12
% EndTime: 2022-01-20 09:12:13
% DurationCPUTime: 0.33s
% Computational Cost: add. (350->100), mult. (708->150), div. (0->0), fcn. (636->8), ass. (0->46)
t38 = cos(pkin(9));
t54 = t38 ^ 2;
t37 = sin(pkin(7));
t28 = pkin(1) * t37 + qJ(3);
t53 = 0.2e1 * t28;
t35 = sin(pkin(9));
t41 = sin(qJ(5));
t42 = cos(qJ(5));
t20 = t35 * t42 + t38 * t41;
t36 = sin(pkin(8));
t12 = t20 * t36;
t19 = -t35 * t41 + t38 * t42;
t13 = t19 * t36;
t46 = -Ifges(6,5) * t13 + Ifges(6,6) * t12;
t16 = (pkin(4) * t35 + t28) * t36;
t47 = t38 * t36;
t48 = t35 * t36;
t17 = mrSges(5,1) * t48 + mrSges(5,2) * t47;
t4 = t12 * mrSges(6,1) + t13 * mrSges(6,2);
t52 = -m(6) * t16 - t17 - t4;
t51 = t13 ^ 2;
t50 = m(5) + m(6);
t39 = cos(pkin(8));
t49 = t28 * t39;
t40 = cos(pkin(7));
t29 = -pkin(1) * t40 - pkin(2);
t18 = -pkin(3) * t39 - qJ(4) * t36 + t29;
t7 = t35 * t18 + t38 * t49;
t45 = t35 ^ 2 + t54;
t32 = t36 ^ 2;
t34 = t39 ^ 2;
t44 = t32 + t34;
t30 = t36 * mrSges(4,2);
t27 = t28 ^ 2;
t24 = t32 * t27;
t22 = -mrSges(5,1) * t39 - mrSges(5,3) * t47;
t21 = mrSges(5,2) * t39 - mrSges(5,3) * t48;
t15 = t38 * t18;
t9 = -mrSges(6,1) * t39 - mrSges(6,3) * t13;
t8 = mrSges(6,2) * t39 - mrSges(6,3) * t12;
t6 = -t35 * t49 + t15;
t5 = -pkin(6) * t48 + t7;
t3 = -pkin(6) * t47 + t15 + (-t28 * t35 - pkin(4)) * t39;
t2 = t3 * t41 + t42 * t5;
t1 = t3 * t42 - t41 * t5;
t10 = [Ifges(6,1) * t51 + 0.2e1 * t1 * t9 + 0.2e1 * t16 * t4 + 0.2e1 * t2 * t8 + 0.2e1 * t7 * t21 + 0.2e1 * t6 * t22 + 0.2e1 * t29 * t30 + Ifges(2,3) + Ifges(3,3) - (0.2e1 * Ifges(6,4) * t13 - Ifges(6,2) * t12) * t12 + 0.2e1 * (t40 * mrSges(3,1) - t37 * mrSges(3,2)) * pkin(1) + t44 * mrSges(4,3) * t53 + (-0.2e1 * t29 * mrSges(4,1) + (Ifges(6,3) + Ifges(4,2) + Ifges(5,3)) * t39 + 0.2e1 * t46) * t39 + m(6) * (t1 ^ 2 + t16 ^ 2 + t2 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2 + t24) + m(4) * (t27 * t34 + t29 ^ 2 + t24) + m(3) * (t37 ^ 2 + t40 ^ 2) * pkin(1) ^ 2 + (t17 * t53 + (Ifges(5,1) * t54 + Ifges(4,1) + (-0.2e1 * Ifges(5,4) * t38 + Ifges(5,2) * t35) * t35) * t36 + 0.2e1 * (-Ifges(5,5) * t38 + Ifges(5,6) * t35 + Ifges(4,4)) * t39) * t36; m(6) * (-t1 * t12 + t13 * t2) + t13 * t8 - t12 * t9 + (t38 * t21 - t35 * t22 + m(5) * (-t35 * t6 + t38 * t7 - t49)) * t36 + t52 * t39; m(3) + m(4) * t44 + m(5) * (t32 * t45 + t34) + m(6) * (t12 ^ 2 + t34 + t51); -t39 * mrSges(4,1) + t19 * t9 + t20 * t8 + t35 * t21 + t38 * t22 + t30 + m(6) * (t1 * t19 + t2 * t20) + m(5) * (t35 * t7 + t38 * t6) + m(4) * t29; m(6) * (-t12 * t19 + t13 * t20); m(4) + m(5) * t45 + m(6) * (t19 ^ 2 + t20 ^ 2); m(5) * t36 * t28 - t52; -t50 * t39; 0; t50; mrSges(6,1) * t1 - mrSges(6,2) * t2 - Ifges(6,3) * t39 - t46; -t4; mrSges(6,1) * t19 - mrSges(6,2) * t20; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;

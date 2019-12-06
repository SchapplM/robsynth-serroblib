% Calculate joint inertia matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:48
% EndTime: 2019-12-05 15:21:49
% DurationCPUTime: 0.29s
% Computational Cost: add. (275->94), mult. (618->139), div. (0->0), fcn. (561->6), ass. (0->42)
t34 = cos(pkin(9));
t49 = t34 ^ 2;
t48 = 2 * qJ(3);
t32 = sin(pkin(9));
t36 = sin(qJ(5));
t37 = cos(qJ(5));
t18 = t37 * t32 + t36 * t34;
t33 = sin(pkin(8));
t12 = t18 * t33;
t17 = -t36 * t32 + t37 * t34;
t13 = t17 * t33;
t42 = -Ifges(6,5) * t13 + Ifges(6,6) * t12;
t43 = t34 * t33;
t44 = t32 * t33;
t14 = mrSges(5,1) * t44 + mrSges(5,2) * t43;
t21 = (pkin(4) * t32 + qJ(3)) * t33;
t3 = t12 * mrSges(6,1) + t13 * mrSges(6,2);
t47 = -m(6) * t21 - t14 - t3;
t46 = t13 ^ 2;
t45 = m(5) + m(6);
t35 = cos(pkin(8));
t22 = -t35 * pkin(3) - t33 * qJ(4) - pkin(2);
t39 = qJ(3) * t35;
t9 = t32 * t22 + t34 * t39;
t41 = t32 ^ 2 + t49;
t29 = t33 ^ 2;
t31 = t35 ^ 2;
t40 = t29 + t31;
t38 = qJ(3) ^ 2;
t27 = t33 * mrSges(4,2);
t26 = t29 * t38;
t20 = -t35 * mrSges(5,1) - mrSges(5,3) * t43;
t19 = t35 * mrSges(5,2) - mrSges(5,3) * t44;
t16 = t34 * t22;
t8 = -t32 * t39 + t16;
t7 = -t35 * mrSges(6,1) - t13 * mrSges(6,3);
t6 = t35 * mrSges(6,2) - t12 * mrSges(6,3);
t5 = -pkin(6) * t44 + t9;
t4 = -pkin(6) * t43 + t16 + (-qJ(3) * t32 - pkin(4)) * t35;
t2 = t36 * t4 + t37 * t5;
t1 = -t36 * t5 + t37 * t4;
t10 = [m(2) + m(3) + m(4) * t40 + m(5) * (t41 * t29 + t31) + m(6) * (t12 ^ 2 + t31 + t46); m(6) * (-t1 * t12 + t2 * t13) + t13 * t6 - t12 * t7 + (t34 * t19 - t32 * t20 + m(5) * (-t32 * t8 + t34 * t9 - t39)) * t33 + t47 * t35; Ifges(6,1) * t46 - 0.2e1 * pkin(2) * t27 + 0.2e1 * t1 * t7 + 0.2e1 * t9 * t19 + 0.2e1 * t2 * t6 + 0.2e1 * t8 * t20 + 0.2e1 * t21 * t3 + Ifges(3,3) - (0.2e1 * Ifges(6,4) * t13 - Ifges(6,2) * t12) * t12 + t40 * mrSges(4,3) * t48 + m(6) * (t1 ^ 2 + t2 ^ 2 + t21 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2 + t26) + m(4) * (pkin(2) ^ 2 + t31 * t38 + t26) + (0.2e1 * pkin(2) * mrSges(4,1) + (Ifges(6,3) + Ifges(5,3) + Ifges(4,2)) * t35 + 0.2e1 * t42) * t35 + (t14 * t48 + (Ifges(5,1) * t49 + Ifges(4,1) + (-0.2e1 * Ifges(5,4) * t34 + Ifges(5,2) * t32) * t32) * t33 + 0.2e1 * (-Ifges(5,5) * t34 + Ifges(5,6) * t32 + Ifges(4,4)) * t35) * t33; m(6) * (-t17 * t12 + t18 * t13); -m(4) * pkin(2) - t35 * mrSges(4,1) + t17 * t7 + t18 * t6 + t32 * t19 + t34 * t20 + t27 + m(6) * (t17 * t1 + t18 * t2) + m(5) * (t32 * t9 + t34 * t8); m(4) + m(5) * t41 + m(6) * (t17 ^ 2 + t18 ^ 2); -t45 * t35; m(5) * t33 * qJ(3) - t47; 0; t45; -t3; t1 * mrSges(6,1) - t2 * mrSges(6,2) - Ifges(6,3) * t35 - t42; t17 * mrSges(6,1) - t18 * mrSges(6,2); 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;

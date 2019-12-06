% Calculate joint inertia matrix for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:06
% EndTime: 2019-12-05 17:06:07
% DurationCPUTime: 0.22s
% Computational Cost: add. (221->70), mult. (430->98), div. (0->0), fcn. (256->6), ass. (0->45)
t32 = cos(qJ(5));
t28 = t32 ^ 2;
t29 = sin(qJ(5));
t27 = t29 ^ 2;
t49 = mrSges(6,3) * t27;
t22 = pkin(8) * t49;
t48 = mrSges(6,3) * t28;
t23 = pkin(8) * t48;
t15 = -t32 * mrSges(6,1) + t29 * mrSges(6,2);
t52 = pkin(4) * t15;
t55 = t22 + t23 - t52;
t33 = cos(qJ(4));
t51 = t33 * pkin(3);
t20 = -pkin(4) - t51;
t10 = t20 * t15;
t30 = sin(qJ(4));
t19 = t30 * pkin(3) + pkin(8);
t13 = t19 * t49;
t14 = t19 * t48;
t24 = mrSges(5,1) * t51;
t54 = t10 + t13 + t14 + t24;
t31 = sin(qJ(3));
t53 = pkin(2) * t31;
t34 = cos(qJ(3));
t21 = t34 * pkin(2) + pkin(3);
t9 = t30 * t21 + t33 * t53;
t50 = t9 * mrSges(5,2);
t47 = t30 * mrSges(5,2);
t46 = Ifges(6,5) * t29 + Ifges(6,6) * t32;
t45 = t27 + t28;
t44 = pkin(3) * t47;
t43 = Ifges(6,2) * t28 + Ifges(5,3) + (Ifges(6,1) * t29 + 0.2e1 * Ifges(6,4) * t32) * t29;
t42 = t45 * t19;
t41 = Ifges(4,3) + t43;
t40 = -mrSges(6,1) * t29 - mrSges(6,2) * t32;
t8 = t33 * t21 - t30 * t53;
t39 = (t34 * mrSges(4,1) - t31 * mrSges(4,2)) * pkin(2);
t6 = -pkin(4) - t8;
t1 = t6 * t15;
t7 = pkin(8) + t9;
t2 = t7 * t49;
t3 = t7 * t48;
t4 = t8 * mrSges(5,1);
t38 = t1 + t2 + t3 + t4 + t43 - t50;
t5 = [m(6) * t45 + m(2) + m(3) + m(4) + m(5); 0; -0.2e1 * t50 + Ifges(3,3) + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 + 0.2e1 * t4 + 0.2e1 * t39 + m(6) * (t45 * t7 ^ 2 + t6 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t31 ^ 2 + t34 ^ 2) * pkin(2) ^ 2 + t41; 0; t38 + m(6) * (t20 * t6 + t42 * t7) + (-t47 + m(5) * (t30 * t9 + t33 * t8)) * pkin(3) + t39 + Ifges(4,3) + t54; -0.2e1 * t44 + 0.2e1 * t10 + 0.2e1 * t13 + 0.2e1 * t14 + 0.2e1 * t24 + m(6) * (t45 * t19 ^ 2 + t20 ^ 2) + m(5) * (t30 ^ 2 + t33 ^ 2) * pkin(3) ^ 2 + t41; 0; m(6) * (t45 * t7 * pkin(8) - pkin(4) * t6) + t38 + t55; -t44 + m(6) * (-pkin(4) * t20 + pkin(8) * t42) + t43 + t54 + t55; 0.2e1 * t23 + 0.2e1 * t22 - 0.2e1 * t52 + m(6) * (t45 * pkin(8) ^ 2 + pkin(4) ^ 2) + t43; -t15; t40 * t7 + t46; t19 * t40 + t46; pkin(8) * t40 + t46; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;

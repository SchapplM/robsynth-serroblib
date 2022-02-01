% Calculate joint inertia matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:48
% EndTime: 2022-01-23 09:18:48
% DurationCPUTime: 0.27s
% Computational Cost: add. (381->91), mult. (704->124), div. (0->0), fcn. (596->8), ass. (0->43)
t35 = sin(pkin(9));
t37 = cos(pkin(9));
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t22 = t41 * t35 + t39 * t37;
t58 = t22 ^ 2;
t57 = t37 ^ 2;
t21 = -t39 * t35 + t41 * t37;
t5 = -t21 * mrSges(6,1) + t22 * mrSges(6,2);
t56 = 0.2e1 * t5;
t36 = sin(pkin(8));
t55 = pkin(1) * t36;
t54 = t37 * pkin(4);
t38 = cos(pkin(8));
t29 = t38 * pkin(1) + pkin(2);
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t14 = t42 * t29 - t40 * t55;
t53 = t14 * mrSges(4,1);
t15 = t40 * t29 + t42 * t55;
t52 = t15 * mrSges(4,2);
t51 = Ifges(6,5) * t22 + Ifges(6,6) * t21;
t50 = t35 ^ 2 + t57;
t49 = 2 * mrSges(6,3);
t24 = -t37 * mrSges(5,1) + t35 * mrSges(5,2);
t48 = t50 * qJ(4);
t13 = -pkin(3) - t14;
t47 = Ifges(6,1) * t58 + Ifges(5,2) * t57 + Ifges(4,3) + (Ifges(5,1) * t35 + 0.2e1 * Ifges(5,4) * t37) * t35 + (0.2e1 * Ifges(6,4) * t22 + Ifges(6,2) * t21) * t21;
t46 = 0.2e1 * t50 * mrSges(5,3);
t45 = t24 + t5;
t32 = t37 * pkin(7);
t30 = -pkin(3) - t54;
t25 = t37 * qJ(4) + t32;
t23 = (-pkin(7) - qJ(4)) * t35;
t12 = qJ(4) + t15;
t10 = t13 - t54;
t9 = t37 * t12 + t32;
t8 = (-pkin(7) - t12) * t35;
t7 = t39 * t23 + t41 * t25;
t6 = t41 * t23 - t39 * t25;
t2 = t39 * t8 + t41 * t9;
t1 = -t39 * t9 + t41 * t8;
t3 = [0.2e1 * t53 - 0.2e1 * t52 + t10 * t56 + 0.2e1 * t13 * t24 + Ifges(2,3) + Ifges(3,3) + (-t1 * t22 + t2 * t21) * t49 + t12 * t46 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(5) * (t50 * t12 ^ 2 + t13 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2) + t47 + (0.2e1 * t38 * mrSges(3,1) - 0.2e1 * t36 * mrSges(3,2) + m(3) * (t36 ^ 2 + t38 ^ 2) * pkin(1)) * pkin(1); m(6) * (t1 * t21 + t2 * t22); m(3) + m(4) + m(5) * t50 + m(6) * (t21 ^ 2 + t58); t53 - t52 + (t10 + t30) * t5 + (t13 - pkin(3)) * t24 + m(6) * (t6 * t1 + t30 * t10 + t7 * t2) + m(5) * (-pkin(3) * t13 + t12 * t48) + ((-t1 - t6) * t22 + (t2 + t7) * t21) * mrSges(6,3) + (t50 * t12 + t48) * mrSges(5,3) + t47; m(6) * (t6 * t21 + t7 * t22); -0.2e1 * pkin(3) * t24 + t30 * t56 + (t7 * t21 - t6 * t22) * t49 + qJ(4) * t46 + m(6) * (t30 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t50 * qJ(4) ^ 2 + pkin(3) ^ 2) + t47; m(5) * t13 + m(6) * t10 + t45; 0; -m(5) * pkin(3) + m(6) * t30 + t45; m(5) + m(6); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t51; -t5; t6 * mrSges(6,1) - t7 * mrSges(6,2) + t51; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;

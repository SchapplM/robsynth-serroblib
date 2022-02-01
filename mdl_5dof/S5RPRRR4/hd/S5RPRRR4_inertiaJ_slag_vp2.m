% Calculate joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:27
% EndTime: 2022-01-23 09:34:28
% DurationCPUTime: 0.23s
% Computational Cost: add. (357->78), mult. (642->109), div. (0->0), fcn. (458->8), ass. (0->51)
t38 = cos(qJ(5));
t32 = t38 ^ 2;
t35 = sin(qJ(5));
t31 = t35 ^ 2;
t56 = mrSges(6,3) * t31;
t26 = pkin(8) * t56;
t55 = mrSges(6,3) * t32;
t27 = pkin(8) * t55;
t20 = -t38 * mrSges(6,1) + t35 * mrSges(6,2);
t59 = pkin(4) * t20;
t62 = t26 + t27 - t59;
t39 = cos(qJ(4));
t58 = t39 * pkin(3);
t25 = -pkin(4) - t58;
t15 = t25 * t20;
t36 = sin(qJ(4));
t24 = t36 * pkin(3) + pkin(8);
t18 = t24 * t56;
t19 = t24 * t55;
t28 = mrSges(5,1) * t58;
t61 = t15 + t18 + t19 + t28;
t33 = sin(pkin(9));
t60 = pkin(1) * t33;
t34 = cos(pkin(9));
t23 = t34 * pkin(1) + pkin(2);
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t13 = t40 * t23 - t37 * t60;
t12 = pkin(3) + t13;
t14 = t37 * t23 + t40 * t60;
t9 = t36 * t12 + t39 * t14;
t57 = t9 * mrSges(5,2);
t54 = t13 * mrSges(4,1);
t53 = t14 * mrSges(4,2);
t52 = t36 * mrSges(5,2);
t51 = Ifges(6,5) * t35 + Ifges(6,6) * t38;
t50 = t31 + t32;
t49 = pkin(3) * t52;
t48 = Ifges(6,2) * t32 + Ifges(5,3) + (Ifges(6,1) * t35 + 0.2e1 * Ifges(6,4) * t38) * t35;
t47 = t50 * t24;
t46 = Ifges(4,3) + t48;
t45 = -mrSges(6,1) * t35 - mrSges(6,2) * t38;
t8 = t39 * t12 - t36 * t14;
t6 = -pkin(4) - t8;
t1 = t6 * t20;
t7 = pkin(8) + t9;
t2 = t7 * t56;
t3 = t7 * t55;
t4 = t8 * mrSges(5,1);
t44 = t1 + t2 + t3 + t4 + t48 - t57;
t5 = [0.2e1 * t54 - 0.2e1 * t53 - 0.2e1 * t57 + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 + 0.2e1 * t4 + m(6) * (t50 * t7 ^ 2 + t6 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t13 ^ 2 + t14 ^ 2) + t46 + (0.2e1 * t34 * mrSges(3,1) - 0.2e1 * t33 * mrSges(3,2) + m(3) * (t33 ^ 2 + t34 ^ 2) * pkin(1)) * pkin(1); 0; m(6) * t50 + m(3) + m(4) + m(5); t44 + m(6) * (t25 * t6 + t7 * t47) + (m(5) * (t36 * t9 + t39 * t8) - t52) * pkin(3) + Ifges(4,3) + t54 - t53 + t61; 0; -0.2e1 * t49 + 0.2e1 * t15 + 0.2e1 * t18 + 0.2e1 * t19 + 0.2e1 * t28 + m(6) * (t50 * t24 ^ 2 + t25 ^ 2) + m(5) * (t36 ^ 2 + t39 ^ 2) * pkin(3) ^ 2 + t46; m(6) * (t50 * t7 * pkin(8) - pkin(4) * t6) + t44 + t62; 0; -t49 + m(6) * (-pkin(4) * t25 + pkin(8) * t47) + t48 + t61 + t62; 0.2e1 * t27 + 0.2e1 * t26 - 0.2e1 * t59 + m(6) * (t50 * pkin(8) ^ 2 + pkin(4) ^ 2) + t48; t45 * t7 + t51; -t20; t45 * t24 + t51; t45 * pkin(8) + t51; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;

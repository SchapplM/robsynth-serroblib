% Calculate kinetic energy for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:36
% EndTime: 2022-01-20 12:01:37
% DurationCPUTime: 1.07s
% Computational Cost: add. (2125->132), mult. (3094->195), div. (0->0), fcn. (2460->10), ass. (0->51)
t49 = V_base(5) * pkin(5) + V_base(1);
t50 = -V_base(4) * pkin(5) + V_base(2);
t57 = sin(qJ(1));
t62 = cos(qJ(1));
t40 = -t49 * t57 + t62 * t50;
t44 = t57 * V_base(5) + t62 * V_base(4);
t52 = V_base(6) + qJD(1);
t34 = pkin(1) * t52 - pkin(6) * t44 + t40;
t41 = t62 * t49 + t57 * t50;
t43 = -t57 * V_base(4) + t62 * V_base(5);
t37 = pkin(6) * t43 + t41;
t56 = sin(qJ(2));
t61 = cos(qJ(2));
t25 = t61 * t34 - t37 * t56;
t39 = t43 * t56 + t44 * t61;
t51 = qJD(2) + t52;
t19 = pkin(2) * t51 - pkin(7) * t39 + t25;
t26 = t56 * t34 + t61 * t37;
t38 = t43 * t61 - t44 * t56;
t22 = pkin(7) * t38 + t26;
t55 = sin(qJ(3));
t60 = cos(qJ(3));
t12 = t55 * t19 + t60 * t22;
t48 = qJD(3) + t51;
t10 = pkin(8) * t48 + t12;
t29 = t38 * t60 - t55 * t39;
t30 = t38 * t55 + t39 * t60;
t42 = -pkin(1) * t43 + V_base(3);
t31 = -pkin(2) * t38 + t42;
t15 = -pkin(3) * t29 - pkin(8) * t30 + t31;
t54 = sin(qJ(4));
t59 = cos(qJ(4));
t6 = t59 * t10 + t54 * t15;
t5 = -t10 * t54 + t59 * t15;
t11 = t19 * t60 - t55 * t22;
t28 = qJD(4) - t29;
t9 = -pkin(3) * t48 - t11;
t63 = V_base(3) ^ 2;
t58 = cos(qJ(5));
t53 = sin(qJ(5));
t27 = qJD(5) + t28;
t24 = t30 * t59 + t48 * t54;
t23 = -t30 * t54 + t48 * t59;
t17 = t23 * t53 + t24 * t58;
t16 = t23 * t58 - t24 * t53;
t7 = -pkin(4) * t23 + t9;
t4 = pkin(9) * t23 + t6;
t3 = pkin(4) * t28 - pkin(9) * t24 + t5;
t2 = t3 * t53 + t4 * t58;
t1 = t3 * t58 - t4 * t53;
t8 = m(2) * (t40 ^ 2 + t41 ^ 2 + t63) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t63) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t42 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t40 * mrSges(2,1) - t41 * mrSges(2,2) + Ifges(2,3) * t52 / 0.2e1) * t52 + (t25 * mrSges(3,1) - t26 * mrSges(3,2) + Ifges(3,3) * t51 / 0.2e1) * t51 + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t40 * mrSges(2,3) + Ifges(2,5) * t52 + Ifges(2,1) * t44 / 0.2e1) * t44 + (t42 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,5) * t51 + Ifges(3,1) * t39 / 0.2e1) * t39 + (t31 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t30 / 0.2e1) * t30 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t28 + Ifges(5,1) * t24 / 0.2e1) * t24 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t27 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t41 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,6) * t52 + Ifges(2,2) * t43 / 0.2e1) * t43 + (-t42 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t39 + Ifges(3,6) * t51 + Ifges(3,2) * t38 / 0.2e1) * t38 + (-t31 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,6) * t48 + Ifges(4,2) * t29 / 0.2e1) * t29 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,6) * t28 + Ifges(5,2) * t23 / 0.2e1) * t23 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t27 + Ifges(6,2) * t16 / 0.2e1) * t16;
T = t8;

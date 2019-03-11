% Calculate kinetic energy for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:48
% EndTime: 2019-03-09 20:53:49
% DurationCPUTime: 0.82s
% Computational Cost: add. (2265->149), mult. (2885->196), div. (0->0), fcn. (2240->8), ass. (0->50)
t63 = cos(qJ(4));
t62 = pkin(4) + qJ(6);
t55 = sin(qJ(1));
t58 = cos(qJ(1));
t43 = -t55 * V_base(4) + t58 * V_base(5);
t44 = t55 * V_base(5) + t58 * V_base(4);
t32 = -pkin(1) * t43 - pkin(7) * t44 + V_base(3);
t48 = V_base(5) * pkin(6) + V_base(1);
t49 = -V_base(4) * pkin(6) + V_base(2);
t39 = t58 * t48 + t55 * t49;
t51 = V_base(6) + qJD(1);
t35 = pkin(7) * t51 + t39;
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t24 = t57 * t32 - t35 * t54;
t37 = t44 * t57 + t51 * t54;
t42 = qJD(2) - t43;
t18 = pkin(2) * t42 - pkin(8) * t37 + t24;
t25 = t54 * t32 + t57 * t35;
t36 = -t44 * t54 + t51 * t57;
t21 = pkin(8) * t36 + t25;
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t14 = t53 * t18 + t56 * t21;
t41 = qJD(3) + t42;
t12 = pkin(9) * t41 + t14;
t27 = t36 * t56 - t37 * t53;
t28 = t36 * t53 + t37 * t56;
t38 = -t55 * t48 + t49 * t58;
t34 = -pkin(1) * t51 - t38;
t29 = -pkin(2) * t36 + t34;
t16 = -pkin(3) * t27 - pkin(9) * t28 + t29;
t52 = sin(qJ(4));
t7 = t63 * t12 + t52 * t16;
t13 = t18 * t56 - t53 * t21;
t26 = qJD(4) - t27;
t5 = -qJ(5) * t26 - t7;
t6 = -t52 * t12 + t63 * t16;
t11 = -pkin(3) * t41 - t13;
t61 = qJD(5) - t6;
t23 = t63 * t28 + t52 * t41;
t60 = -qJ(5) * t23 + t11;
t59 = V_base(3) ^ 2;
t22 = t28 * t52 - t63 * t41;
t8 = pkin(4) * t22 + t60;
t4 = -t26 * pkin(4) + t61;
t3 = t62 * t22 + t60;
t2 = -pkin(5) * t22 + qJD(6) - t5;
t1 = t23 * pkin(5) - t62 * t26 + t61;
t9 = m(2) * (t38 ^ 2 + t39 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t34 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t38 * mrSges(2,1) - t39 * mrSges(2,2) + Ifges(2,3) * t51 / 0.2e1) * t51 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t42 / 0.2e1) * t42 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t38 * mrSges(2,3) + Ifges(2,5) * t51 + Ifges(2,1) * t44 / 0.2e1) * t44 + (t34 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,5) * t42 + Ifges(3,1) * t37 / 0.2e1) * t37 + (t29 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t28 / 0.2e1) * t28 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t39 * mrSges(2,3) + Ifges(2,4) * t44 + Ifges(2,6) * t51 + Ifges(2,2) * t43 / 0.2e1) * t43 + (-t34 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,6) * t42 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t29 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t28 + Ifges(4,6) * t41 + Ifges(4,2) * t27 / 0.2e1) * t27 + (t6 * mrSges(5,1) - t7 * mrSges(5,2) + t4 * mrSges(6,2) + t2 * mrSges(7,2) - t5 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t26) * t26 + (t4 * mrSges(6,1) + t1 * mrSges(7,1) + t11 * mrSges(5,2) - t3 * mrSges(7,2) - t6 * mrSges(5,3) - t8 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t23 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t26) * t23 + (t11 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(7,1) - t8 * mrSges(6,2) - t7 * mrSges(5,3) + t3 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t22 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t26 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t23) * t22;
T  = t9;

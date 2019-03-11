% Calculate kinetic energy for
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:18
% EndTime: 2019-03-09 15:28:19
% DurationCPUTime: 0.85s
% Computational Cost: add. (1847->151), mult. (2353->196), div. (0->0), fcn. (1776->8), ass. (0->51)
t67 = -pkin(4) - pkin(9);
t66 = cos(qJ(1));
t58 = sin(qJ(1));
t44 = -t58 * V_base(4) + t66 * V_base(5);
t45 = t58 * V_base(5) + t66 * V_base(4);
t31 = -pkin(1) * t44 - pkin(7) * t45 + V_base(3);
t50 = V_base(5) * pkin(6) + V_base(1);
t51 = -V_base(4) * pkin(6) + V_base(2);
t40 = t66 * t50 + t58 * t51;
t54 = V_base(6) + qJD(1);
t35 = pkin(7) * t54 + t40;
t57 = sin(qJ(2));
t61 = cos(qJ(2));
t21 = t61 * t31 - t35 * t57;
t38 = t45 * t61 + t54 * t57;
t43 = qJD(2) - t44;
t15 = pkin(2) * t43 - pkin(8) * t38 + t21;
t22 = t57 * t31 + t61 * t35;
t37 = -t45 * t57 + t54 * t61;
t18 = pkin(8) * t37 + t22;
t56 = sin(qJ(3));
t60 = cos(qJ(3));
t12 = t56 * t15 + t60 * t18;
t39 = -t58 * t50 + t66 * t51;
t42 = qJD(3) + t43;
t10 = t42 * qJ(4) + t12;
t34 = -t54 * pkin(1) - t39;
t11 = t15 * t60 - t56 * t18;
t28 = -t37 * pkin(2) + t34;
t65 = qJD(4) - t11;
t27 = t37 * t56 + t38 * t60;
t26 = -t60 * t37 + t38 * t56;
t7 = -qJ(5) * t26 - t10;
t13 = t26 * pkin(3) - t27 * qJ(4) + t28;
t64 = qJD(5) - t13;
t63 = -qJ(5) * t27 + t65;
t62 = V_base(3) ^ 2;
t59 = cos(qJ(6));
t55 = sin(qJ(6));
t25 = qJD(6) + t27;
t20 = t26 * t59 - t42 * t55;
t19 = -t26 * t55 - t42 * t59;
t9 = -pkin(3) * t42 + t65;
t8 = -pkin(4) * t26 + t64;
t6 = pkin(5) * t42 - t7;
t5 = (-pkin(3) - pkin(4)) * t42 + t63;
t4 = (-pkin(3) + t67) * t42 + t63;
t3 = pkin(5) * t27 + t67 * t26 + t64;
t2 = t3 * t55 + t4 * t59;
t1 = t3 * t59 - t4 * t55;
t14 = m(2) * (t39 ^ 2 + t40 ^ 2 + t62) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t62) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t13 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t39 * mrSges(2,1) - t40 * mrSges(2,2) + Ifges(2,3) * t54 / 0.2e1) * t54 + (t21 * mrSges(3,1) - t22 * mrSges(3,2) + Ifges(3,3) * t43 / 0.2e1) * t43 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t39 * mrSges(2,3) + Ifges(2,5) * t54 + Ifges(2,1) * t45 / 0.2e1) * t45 + (t34 * mrSges(3,2) - t21 * mrSges(3,3) + Ifges(3,5) * t43 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t6 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t40 * mrSges(2,3) + Ifges(2,4) * t45 + Ifges(2,6) * t54 + Ifges(2,2) * t44 / 0.2e1) * t44 + (-t34 * mrSges(3,1) + t22 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t43 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t6 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t25 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t11 * mrSges(4,1) - t9 * mrSges(5,1) - t7 * mrSges(6,1) - t12 * mrSges(4,2) + t5 * mrSges(6,2) + t10 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t42) * t42 + (t8 * mrSges(6,1) + t28 * mrSges(4,2) + t9 * mrSges(5,2) - t11 * mrSges(4,3) - t13 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t27 + (Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * t42) * t27 + (t28 * mrSges(4,1) + t13 * mrSges(5,1) - t10 * mrSges(5,2) + t8 * mrSges(6,2) - t12 * mrSges(4,3) - t7 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t26 + (-Ifges(6,5) - Ifges(4,6) + Ifges(5,6)) * t42 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t27) * t26;
T  = t14;

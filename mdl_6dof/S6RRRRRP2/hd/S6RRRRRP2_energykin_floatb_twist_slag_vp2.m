% Calculate kinetic energy for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:35
% EndTime: 2019-03-10 01:00:36
% DurationCPUTime: 1.13s
% Computational Cost: add. (3385->152), mult. (4419->213), div. (0->0), fcn. (3572->10), ass. (0->54)
t63 = sin(qJ(1));
t67 = cos(qJ(1));
t51 = t63 * V_base(5) + t67 * V_base(4);
t58 = V_base(6) + qJD(1);
t62 = sin(qJ(2));
t66 = cos(qJ(2));
t42 = -t51 * t62 + t58 * t66;
t43 = t51 * t66 + t58 * t62;
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t33 = t42 * t65 - t43 * t61;
t34 = t42 * t61 + t43 * t65;
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t23 = t33 * t64 - t34 * t60;
t24 = t33 * t60 + t34 * t64;
t55 = V_base(5) * pkin(6) + V_base(1);
t56 = -V_base(4) * pkin(6) + V_base(2);
t44 = -t63 * t55 + t56 * t67;
t40 = -pkin(1) * t58 - t44;
t35 = -pkin(2) * t42 + t40;
t25 = -pkin(3) * t33 + t35;
t12 = -pkin(4) * t23 - pkin(10) * t24 + t25;
t59 = sin(qJ(5));
t69 = cos(qJ(5));
t50 = -t63 * V_base(4) + t67 * V_base(5);
t38 = -pkin(1) * t50 - pkin(7) * t51 + V_base(3);
t45 = t67 * t55 + t63 * t56;
t41 = pkin(7) * t58 + t45;
t31 = t66 * t38 - t41 * t62;
t49 = qJD(2) - t50;
t28 = pkin(2) * t49 - pkin(8) * t43 + t31;
t32 = t62 * t38 + t66 * t41;
t30 = pkin(8) * t42 + t32;
t18 = t65 * t28 - t30 * t61;
t48 = qJD(3) + t49;
t14 = pkin(3) * t48 - pkin(9) * t34 + t18;
t19 = t61 * t28 + t65 * t30;
t17 = pkin(9) * t33 + t19;
t10 = t60 * t14 + t64 * t17;
t47 = qJD(4) + t48;
t8 = pkin(10) * t47 + t10;
t4 = t59 * t12 + t69 * t8;
t9 = t14 * t64 - t60 * t17;
t7 = -pkin(4) * t47 - t9;
t3 = t12 * t69 - t59 * t8;
t68 = V_base(3) ^ 2;
t22 = qJD(5) - t23;
t21 = t24 * t69 + t59 * t47;
t20 = t24 * t59 - t47 * t69;
t5 = pkin(5) * t20 - qJ(6) * t21 + t7;
t2 = qJ(6) * t22 + t4;
t1 = -t22 * pkin(5) + qJD(6) - t3;
t6 = (t31 * mrSges(3,1) - t32 * mrSges(3,2) + Ifges(3,3) * t49 / 0.2e1) * t49 + (t25 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t47 + Ifges(5,1) * t24 / 0.2e1) * t24 + (t18 * mrSges(4,1) - t19 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (t35 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(2) * (t44 ^ 2 + t45 ^ 2 + t68) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t68) / 0.2e1 + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,5) * t58 + Ifges(2,1) * t51 / 0.2e1) * t51 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t22) * t22 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(3) * (t31 ^ 2 + t32 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t35 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t25 ^ 2 + t9 ^ 2) / 0.2e1 + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (t40 * mrSges(3,2) - t31 * mrSges(3,3) + Ifges(3,5) * t49 + Ifges(3,1) * t43 / 0.2e1) * t43 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t25 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,6) * t47 + Ifges(5,2) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,6) * t58 + Ifges(2,2) * t50 / 0.2e1) * t50 + (t7 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t20 + (-Ifges(6,6) + Ifges(7,6)) * t22 + (-Ifges(6,4) + Ifges(7,5)) * t21) * t20 + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + Ifges(2,3) * t58 / 0.2e1) * t58 + (-t40 * mrSges(3,1) + t32 * mrSges(3,3) + Ifges(3,4) * t43 + Ifges(3,6) * t49 + Ifges(3,2) * t42 / 0.2e1) * t42 + (t7 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t21 + (Ifges(7,4) + Ifges(6,5)) * t22) * t21 + (-t35 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t34 + Ifges(4,6) * t48 + Ifges(4,2) * t33 / 0.2e1) * t33;
T  = t6;

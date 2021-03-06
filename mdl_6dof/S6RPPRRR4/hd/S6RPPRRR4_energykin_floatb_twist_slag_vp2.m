% Calculate kinetic energy for
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:27
% EndTime: 2019-03-09 02:25:28
% DurationCPUTime: 0.96s
% Computational Cost: add. (2243->153), mult. (3039->211), div. (0->0), fcn. (2276->10), ass. (0->54)
t64 = sin(qJ(1));
t70 = cos(qJ(1));
t49 = t64 * V_base(4) - t70 * V_base(5);
t50 = t64 * V_base(5) + t70 * V_base(4);
t42 = t49 * pkin(1) - t50 * qJ(2) + V_base(3);
t31 = -pkin(2) * t49 + qJD(3) - t42;
t59 = sin(pkin(10));
t60 = cos(pkin(10));
t40 = t49 * t60 - t50 * t59;
t41 = t49 * t59 + t50 * t60;
t19 = -pkin(3) * t40 - pkin(7) * t41 + t31;
t58 = V_base(6) + qJD(1);
t54 = V_base(5) * pkin(6) + V_base(1);
t55 = -V_base(4) * pkin(6) + V_base(2);
t45 = -t64 * t54 + t70 * t55;
t69 = qJD(2) - t45;
t29 = -t50 * qJ(3) + (-pkin(1) - pkin(2)) * t58 + t69;
t46 = t70 * t54 + t64 * t55;
t44 = t58 * qJ(2) + t46;
t37 = qJ(3) * t49 + t44;
t25 = t59 * t29 + t60 * t37;
t23 = -pkin(7) * t58 + t25;
t63 = sin(qJ(4));
t67 = cos(qJ(4));
t12 = t63 * t19 + t67 * t23;
t39 = qJD(4) - t40;
t10 = pkin(8) * t39 + t12;
t24 = t29 * t60 - t59 * t37;
t22 = pkin(3) * t58 - t24;
t35 = -t63 * t41 - t58 * t67;
t36 = t41 * t67 - t58 * t63;
t15 = -pkin(4) * t35 - pkin(8) * t36 + t22;
t62 = sin(qJ(5));
t66 = cos(qJ(5));
t6 = t66 * t10 + t62 * t15;
t5 = -t10 * t62 + t66 * t15;
t11 = t19 * t67 - t63 * t23;
t32 = qJD(5) - t35;
t9 = -pkin(4) * t39 - t11;
t68 = V_base(3) ^ 2;
t65 = cos(qJ(6));
t61 = sin(qJ(6));
t43 = -t58 * pkin(1) + t69;
t30 = qJD(6) + t32;
t27 = t36 * t66 + t39 * t62;
t26 = -t36 * t62 + t39 * t66;
t17 = t26 * t61 + t27 * t65;
t16 = t26 * t65 - t27 * t61;
t7 = -pkin(5) * t26 + t9;
t4 = pkin(9) * t26 + t6;
t3 = pkin(5) * t32 - pkin(9) * t27 + t5;
t2 = t3 * t61 + t4 * t65;
t1 = t3 * t65 - t4 * t61;
t8 = m(2) * (t45 ^ 2 + t46 ^ 2 + t68) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t68) / 0.2e1 + m(3) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t22 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t31 * mrSges(4,2) - t24 * mrSges(4,3) + Ifges(4,1) * t41 / 0.2e1) * t41 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t39 / 0.2e1) * t39 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t32 / 0.2e1) * t32 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t30 / 0.2e1) * t30 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t31 * mrSges(4,1) + t25 * mrSges(4,3) + Ifges(4,4) * t41 + Ifges(4,2) * t40 / 0.2e1) * t40 + (t22 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t39 + Ifges(5,1) * t36 / 0.2e1) * t36 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t32 + Ifges(6,1) * t27 / 0.2e1) * t27 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t30 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t22 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t36 + Ifges(5,6) * t39 + Ifges(5,2) * t35 / 0.2e1) * t35 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t32 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t30 + Ifges(7,2) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(2,2) + t43 * mrSges(3,2) - t45 * mrSges(2,3) - t42 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t50) * t50 + (V_base(3) * mrSges(2,1) + t42 * mrSges(3,1) - t44 * mrSges(3,2) - t46 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t49 + (-Ifges(2,4) + Ifges(3,5)) * t50) * t49 + (t45 * mrSges(2,1) - t43 * mrSges(3,1) - t24 * mrSges(4,1) - t46 * mrSges(2,2) + t25 * mrSges(4,2) + t44 * mrSges(3,3) - Ifges(4,5) * t41 - Ifges(4,6) * t40 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * t58 + (Ifges(3,4) + Ifges(2,5)) * t50 + (-Ifges(2,6) + Ifges(3,6)) * t49) * t58;
T  = t8;

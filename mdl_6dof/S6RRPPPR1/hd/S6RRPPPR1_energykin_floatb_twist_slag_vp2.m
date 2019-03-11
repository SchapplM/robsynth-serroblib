% Calculate kinetic energy for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:05:58
% EndTime: 2019-03-09 08:05:59
% DurationCPUTime: 1.13s
% Computational Cost: add. (2665->153), mult. (3523->209), div. (0->0), fcn. (2788->10), ass. (0->55)
t71 = -pkin(4) - pkin(5);
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t49 = -t61 * V_base(4) + t64 * V_base(5);
t50 = t61 * V_base(5) + t64 * V_base(4);
t38 = -pkin(1) * t49 - pkin(7) * t50 + V_base(3);
t54 = V_base(5) * pkin(6) + V_base(1);
t55 = -V_base(4) * pkin(6) + V_base(2);
t46 = t64 * t54 + t61 * t55;
t56 = V_base(6) + qJD(1);
t42 = pkin(7) * t56 + t46;
t60 = sin(qJ(2));
t63 = cos(qJ(2));
t29 = t63 * t38 - t42 * t60;
t44 = t50 * t63 + t56 * t60;
t48 = qJD(2) - t49;
t23 = pkin(2) * t48 - qJ(3) * t44 + t29;
t30 = t60 * t38 + t63 * t42;
t43 = -t50 * t60 + t56 * t63;
t26 = qJ(3) * t43 + t30;
t58 = sin(pkin(9));
t70 = cos(pkin(9));
t16 = t58 * t23 + t70 * t26;
t14 = qJ(4) * t48 + t16;
t33 = -t70 * t43 + t44 * t58;
t34 = t58 * t43 + t44 * t70;
t45 = -t61 * t54 + t55 * t64;
t41 = -pkin(1) * t56 - t45;
t35 = -pkin(2) * t43 + qJD(3) + t41;
t18 = pkin(3) * t33 - qJ(4) * t34 + t35;
t57 = sin(pkin(10));
t69 = cos(pkin(10));
t9 = t69 * t14 + t57 * t18;
t15 = t70 * t23 - t58 * t26;
t7 = t33 * qJ(5) + t9;
t8 = -t57 * t14 + t18 * t69;
t68 = pkin(3) * t48 - qJD(4) + t15;
t67 = qJD(5) - t8;
t28 = t34 * t69 + t57 * t48;
t66 = qJ(5) * t28 + t68;
t65 = V_base(3) ^ 2;
t62 = cos(qJ(6));
t59 = sin(qJ(6));
t32 = qJD(6) - t33;
t27 = t34 * t57 - t48 * t69;
t20 = t27 * t59 + t28 * t62;
t19 = t27 * t62 - t28 * t59;
t10 = pkin(4) * t27 - t66;
t6 = -t33 * pkin(4) + t67;
t5 = t27 * t71 + t66;
t4 = pkin(8) * t27 + t7;
t3 = -t28 * pkin(8) + t33 * t71 + t67;
t2 = t3 * t59 + t4 * t62;
t1 = t3 * t62 - t4 * t59;
t11 = m(2) * (t45 ^ 2 + t46 ^ 2 + t65) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t65) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t41 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t45 * mrSges(2,1) - t46 * mrSges(2,2) + Ifges(2,3) * t56 / 0.2e1) * t56 + (t41 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,1) * t44 / 0.2e1) * t44 + (t35 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,1) * t34 / 0.2e1) * t34 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,5) * t56 + Ifges(2,1) * t50 / 0.2e1) * t50 + (-t41 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,2) * t43 / 0.2e1) * t43 + (t5 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t32 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t46 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,6) * t56 + Ifges(2,2) * t49 / 0.2e1) * t49 + (-t5 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t32 + Ifges(7,2) * t19 / 0.2e1) * t19 + (-t68 * mrSges(5,2) + t6 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t28) * t28 + (-t68 * mrSges(5,1) + t10 * mrSges(6,1) - t7 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t27 + (-Ifges(5,4) + Ifges(6,5)) * t28) * t27 + (t29 * mrSges(3,1) + t15 * mrSges(4,1) - t30 * mrSges(3,2) - t16 * mrSges(4,2) + Ifges(3,5) * t44 + Ifges(4,5) * t34 + Ifges(3,6) * t43 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t48) * t48 + (t35 * mrSges(4,1) + t8 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t16 * mrSges(4,3) + t7 * mrSges(6,3) - Ifges(4,4) * t34 - Ifges(4,6) * t48 + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t33 + (Ifges(6,4) + Ifges(5,5)) * t28 + (-Ifges(5,6) + Ifges(6,6)) * t27) * t33;
T  = t11;

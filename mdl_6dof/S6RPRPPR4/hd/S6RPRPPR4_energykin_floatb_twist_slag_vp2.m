% Calculate kinetic energy for
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:40
% EndTime: 2019-03-09 02:46:41
% DurationCPUTime: 1.12s
% Computational Cost: add. (2635->153), mult. (3523->209), div. (0->0), fcn. (2788->10), ass. (0->55)
t72 = -pkin(4) - pkin(5);
t71 = cos(qJ(3));
t63 = sin(qJ(1));
t65 = cos(qJ(1));
t49 = t63 * V_base(4) - t65 * V_base(5);
t50 = t63 * V_base(5) + t65 * V_base(4);
t38 = pkin(1) * t49 - qJ(2) * t50 + V_base(3);
t54 = V_base(5) * pkin(6) + V_base(1);
t55 = -V_base(4) * pkin(6) + V_base(2);
t46 = t65 * t54 + t63 * t55;
t57 = V_base(6) + qJD(1);
t42 = qJ(2) * t57 + t46;
t59 = sin(pkin(9));
t60 = cos(pkin(9));
t29 = t60 * t38 - t42 * t59;
t44 = t50 * t60 + t57 * t59;
t23 = pkin(2) * t49 - pkin(7) * t44 + t29;
t30 = t59 * t38 + t60 * t42;
t43 = -t50 * t59 + t57 * t60;
t26 = pkin(7) * t43 + t30;
t62 = sin(qJ(3));
t16 = t62 * t23 + t71 * t26;
t48 = qJD(3) + t49;
t14 = qJ(4) * t48 + t16;
t33 = -t71 * t43 + t44 * t62;
t34 = t62 * t43 + t44 * t71;
t45 = -t63 * t54 + t55 * t65;
t40 = -pkin(1) * t57 + qJD(2) - t45;
t35 = -pkin(2) * t43 + t40;
t18 = pkin(3) * t33 - qJ(4) * t34 + t35;
t58 = sin(pkin(10));
t70 = cos(pkin(10));
t9 = t70 * t14 + t58 * t18;
t15 = t71 * t23 - t62 * t26;
t7 = t33 * qJ(5) + t9;
t8 = -t58 * t14 + t18 * t70;
t69 = pkin(3) * t48 - qJD(4) + t15;
t68 = qJD(5) - t8;
t28 = t34 * t70 + t58 * t48;
t67 = qJ(5) * t28 + t69;
t66 = V_base(3) ^ 2;
t64 = cos(qJ(6));
t61 = sin(qJ(6));
t32 = qJD(6) - t33;
t27 = t34 * t58 - t48 * t70;
t20 = t27 * t61 + t28 * t64;
t19 = t27 * t64 - t28 * t61;
t10 = pkin(4) * t27 - t67;
t6 = -t33 * pkin(4) + t68;
t5 = t27 * t72 + t67;
t4 = pkin(8) * t27 + t7;
t3 = -t28 * pkin(8) + t33 * t72 + t68;
t2 = t3 * t61 + t4 * t64;
t1 = t3 * t64 - t4 * t61;
t11 = m(2) * (t45 ^ 2 + t46 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t40 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t35 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t69 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t45 * mrSges(2,1) - t46 * mrSges(2,2) + Ifges(2,3) * t57 / 0.2e1) * t57 + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (t40 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,1) * t44 / 0.2e1) * t44 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,5) * t57 + Ifges(2,1) * t50 / 0.2e1) * t50 + (-t40 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,2) * t43 / 0.2e1) * t43 + (t35 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t34 / 0.2e1) * t34 + (t5 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t32 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t5 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t32 + Ifges(7,2) * t19 / 0.2e1) * t19 + (-t69 * mrSges(5,2) + t6 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t28) * t28 + (-t69 * mrSges(5,1) + t10 * mrSges(6,1) - t7 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t27 + (-Ifges(5,4) + Ifges(6,5)) * t28) * t27 + (V_base(3) * mrSges(2,1) + t29 * mrSges(3,1) - t30 * mrSges(3,2) - t46 * mrSges(2,3) - Ifges(2,4) * t50 + Ifges(3,5) * t44 - Ifges(2,6) * t57 + Ifges(3,6) * t43 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t49) * t49 + (t35 * mrSges(4,1) + t8 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t16 * mrSges(4,3) + t7 * mrSges(6,3) - Ifges(4,4) * t34 - Ifges(4,6) * t48 + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t33 + (Ifges(6,4) + Ifges(5,5)) * t28 + (-Ifges(5,6) + Ifges(6,6)) * t27) * t33;
T  = t11;

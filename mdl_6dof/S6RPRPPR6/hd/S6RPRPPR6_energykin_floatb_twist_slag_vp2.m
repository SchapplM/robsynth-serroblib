% Calculate kinetic energy for
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:06
% EndTime: 2019-03-09 02:53:07
% DurationCPUTime: 0.95s
% Computational Cost: add. (2381->153), mult. (3055->209), div. (0->0), fcn. (2284->10), ass. (0->55)
t69 = pkin(1) + pkin(7);
t61 = sin(qJ(1));
t68 = cos(qJ(1));
t47 = t61 * V_base(5) + t68 * V_base(4);
t55 = V_base(6) + qJD(1);
t51 = V_base(5) * pkin(6) + V_base(1);
t52 = -V_base(4) * pkin(6) + V_base(2);
t43 = -t61 * t51 + t52 * t68;
t65 = qJD(2) - t43;
t33 = t47 * pkin(2) - t55 * t69 + t65;
t46 = t61 * V_base(4) - t68 * V_base(5);
t66 = -qJ(2) * t47 + V_base(3);
t35 = t46 * t69 + t66;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t23 = t63 * t33 - t35 * t60;
t42 = t46 * t60 + t55 * t63;
t45 = qJD(3) + t47;
t19 = pkin(3) * t45 - qJ(4) * t42 + t23;
t24 = t60 * t33 + t63 * t35;
t41 = t46 * t63 - t55 * t60;
t22 = qJ(4) * t41 + t24;
t57 = sin(pkin(9));
t67 = cos(pkin(9));
t12 = t57 * t19 + t67 * t22;
t10 = qJ(5) * t45 + t12;
t44 = t68 * t51 + t61 * t52;
t40 = -t55 * qJ(2) - t44;
t36 = -pkin(2) * t46 - t40;
t27 = -pkin(3) * t41 + qJD(4) + t36;
t29 = -t67 * t41 + t42 * t57;
t30 = t57 * t41 + t42 * t67;
t15 = pkin(4) * t29 - qJ(5) * t30 + t27;
t56 = sin(pkin(10));
t58 = cos(pkin(10));
t6 = t58 * t10 + t56 * t15;
t5 = -t10 * t56 + t58 * t15;
t11 = t19 * t67 - t57 * t22;
t9 = -t45 * pkin(4) + qJD(5) - t11;
t64 = V_base(3) ^ 2;
t62 = cos(qJ(6));
t59 = sin(qJ(6));
t38 = -t55 * pkin(1) + t65;
t37 = pkin(1) * t46 + t66;
t28 = qJD(6) + t29;
t26 = t30 * t58 + t45 * t56;
t25 = -t30 * t56 + t45 * t58;
t17 = t25 * t59 + t26 * t62;
t16 = t25 * t62 - t26 * t59;
t7 = -t25 * pkin(5) + t9;
t4 = pkin(8) * t25 + t6;
t3 = pkin(5) * t29 - pkin(8) * t26 + t5;
t2 = t3 * t59 + t4 * t62;
t1 = t3 * t62 - t4 * t59;
t8 = m(2) * (t43 ^ 2 + t44 ^ 2 + t64) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t64) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t36 ^ 2) / 0.2e1 + m(3) * (t37 ^ 2 + t38 ^ 2 + t40 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t27 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t36 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,1) * t42 / 0.2e1) * t42 + (t27 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,1) * t30 / 0.2e1) * t30 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t28 / 0.2e1) * t28 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t36 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t42 + Ifges(4,2) * t41 / 0.2e1) * t41 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,2) * t25 / 0.2e1) * t25 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t28 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t28 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + t38 * mrSges(3,2) - t40 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t55) * t55 + (t38 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) - t37 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t47 + (-Ifges(3,4) + Ifges(2,5)) * t55) * t47 + (t23 * mrSges(4,1) + t11 * mrSges(5,1) - t24 * mrSges(4,2) - t12 * mrSges(5,2) + Ifges(4,5) * t42 + Ifges(5,5) * t30 + Ifges(4,6) * t41 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t45) * t45 + (V_base(3) * mrSges(2,1) + t40 * mrSges(3,1) - t37 * mrSges(3,2) - t44 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t46 + (Ifges(3,5) - Ifges(2,6)) * t55 + (-Ifges(2,4) - Ifges(3,6)) * t47) * t46 + (t27 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t30 + Ifges(6,5) * t26 - Ifges(5,6) * t45 + Ifges(6,6) * t25 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t29) * t29;
T  = t8;

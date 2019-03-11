% Calculate kinetic energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:39
% EndTime: 2019-03-09 15:47:40
% DurationCPUTime: 1.24s
% Computational Cost: add. (4501->158), mult. (6763->221), div. (0->0), fcn. (5696->12), ass. (0->62)
t55 = V_base(5) * pkin(7) + V_base(1);
t56 = -V_base(4) * pkin(7) + V_base(2);
t64 = sin(qJ(1));
t68 = cos(qJ(1));
t48 = -t55 * t64 + t68 * t56;
t57 = V_base(6) + qJD(1);
t60 = cos(pkin(6));
t51 = t64 * V_base(5) + t68 * V_base(4);
t78 = pkin(8) * t51;
t43 = pkin(1) * t57 - t60 * t78 + t48;
t50 = -t64 * V_base(4) + t68 * V_base(5);
t59 = sin(pkin(6));
t46 = -pkin(1) * t50 - t59 * t78 + V_base(3);
t80 = t43 * t60 + t46 * t59;
t49 = t68 * t55 + t64 * t56;
t72 = t50 * t60 + t57 * t59;
t40 = pkin(8) * t72 + t49;
t63 = sin(qJ(2));
t67 = cos(qJ(2));
t30 = -t63 * t40 + t67 * t80;
t41 = -t51 * t63 + t67 * t72;
t79 = pkin(4) + pkin(10);
t32 = -t43 * t59 + t60 * t46;
t42 = t51 * t67 + t63 * t72;
t23 = -pkin(2) * t41 - pkin(9) * t42 + t32;
t31 = t67 * t40 + t80 * t63;
t47 = -t50 * t59 + t57 * t60 + qJD(2);
t26 = pkin(9) * t47 + t31;
t62 = sin(qJ(3));
t66 = cos(qJ(3));
t16 = t66 * t23 - t26 * t62;
t35 = t42 * t66 + t47 * t62;
t39 = qJD(3) - t41;
t12 = pkin(3) * t39 - qJ(4) * t35 + t16;
t17 = t62 * t23 + t66 * t26;
t34 = -t42 * t62 + t47 * t66;
t15 = qJ(4) * t34 + t17;
t58 = sin(pkin(11));
t74 = cos(pkin(11));
t8 = t58 * t12 + t74 * t15;
t6 = -qJ(5) * t39 - t8;
t7 = t12 * t74 - t58 * t15;
t71 = qJD(5) - t7;
t29 = t58 * t34 + t35 * t74;
t25 = -pkin(2) * t47 - t30;
t18 = -pkin(3) * t34 + qJD(4) + t25;
t70 = -qJ(5) * t29 + t18;
t69 = V_base(3) ^ 2;
t65 = cos(qJ(6));
t61 = sin(qJ(6));
t28 = -t34 * t74 + t35 * t58;
t27 = qJD(6) + t29;
t20 = t28 * t61 + t39 * t65;
t19 = t28 * t65 - t39 * t61;
t10 = pkin(4) * t28 + t70;
t9 = t28 * t79 + t70;
t5 = -t39 * pkin(4) + t71;
t4 = -pkin(5) * t28 - t6;
t3 = t29 * pkin(5) - t39 * t79 + t71;
t2 = t3 * t61 + t65 * t9;
t1 = t3 * t65 - t61 * t9;
t11 = m(2) * (t48 ^ 2 + t49 ^ 2 + t69) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t69) / 0.2e1 + m(3) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t25 ^ 2) / 0.2e1 + m(5) * (t18 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t48 * mrSges(2,1) - t49 * mrSges(2,2) + Ifges(2,3) * t57 / 0.2e1) * t57 + (t30 * mrSges(3,1) - t31 * mrSges(3,2) + Ifges(3,3) * t47 / 0.2e1) * t47 + (t25 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,1) * t35 / 0.2e1) * t35 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t48 * mrSges(2,3) + Ifges(2,5) * t57 + Ifges(2,1) * t51 / 0.2e1) * t51 + (t32 * mrSges(3,2) - t30 * mrSges(3,3) + Ifges(3,5) * t47 + Ifges(3,1) * t42 / 0.2e1) * t42 + (-t25 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t35 + Ifges(4,2) * t34 / 0.2e1) * t34 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t27 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t49 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,6) * t57 + Ifges(2,2) * t50 / 0.2e1) * t50 + (-t32 * mrSges(3,1) + t31 * mrSges(3,3) + Ifges(3,4) * t42 + Ifges(3,6) * t47 + Ifges(3,2) * t41 / 0.2e1) * t41 + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t27 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t5 * mrSges(6,1) + t18 * mrSges(5,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t29) * t29 + (t18 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t28 + (-Ifges(5,4) - Ifges(6,6)) * t29) * t28 + (t16 * mrSges(4,1) + t7 * mrSges(5,1) - t17 * mrSges(4,2) - t8 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + Ifges(4,5) * t35 + Ifges(4,6) * t34 + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t39 + (-Ifges(6,4) + Ifges(5,5)) * t29 + (Ifges(6,5) - Ifges(5,6)) * t28) * t39;
T  = t11;

% Calculate kinetic energy for
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR13_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR13_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:50:57
% EndTime: 2019-03-09 23:50:59
% DurationCPUTime: 1.42s
% Computational Cost: add. (4329->158), mult. (6449->223), div. (0->0), fcn. (5406->12), ass. (0->63)
t59 = V_base(5) * pkin(7) + V_base(1);
t60 = -V_base(4) * pkin(7) + V_base(2);
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t52 = -t59 * t68 + t72 * t60;
t61 = V_base(6) + qJD(1);
t63 = cos(pkin(6));
t55 = t68 * V_base(5) + t72 * V_base(4);
t83 = pkin(8) * t55;
t46 = pkin(1) * t61 - t63 * t83 + t52;
t54 = -t68 * V_base(4) + t72 * V_base(5);
t62 = sin(pkin(6));
t49 = -pkin(1) * t54 - t62 * t83 + V_base(3);
t85 = t46 * t63 + t49 * t62;
t53 = t72 * t59 + t68 * t60;
t76 = t54 * t63 + t61 * t62;
t43 = t76 * pkin(8) + t53;
t67 = sin(qJ(2));
t71 = cos(qJ(2));
t28 = -t67 * t43 + t85 * t71;
t44 = -t55 * t67 + t76 * t71;
t84 = -pkin(4) - pkin(5);
t82 = cos(qJ(4));
t32 = -t46 * t62 + t63 * t49;
t45 = t55 * t71 + t76 * t67;
t23 = -pkin(2) * t44 - pkin(9) * t45 + t32;
t29 = t71 * t43 + t85 * t67;
t51 = -t54 * t62 + t61 * t63 + qJD(2);
t27 = pkin(9) * t51 + t29;
t66 = sin(qJ(3));
t70 = cos(qJ(3));
t16 = t66 * t23 + t70 * t27;
t42 = qJD(3) - t44;
t14 = pkin(10) * t42 + t16;
t26 = -pkin(2) * t51 - t28;
t36 = -t45 * t66 + t70 * t51;
t37 = t45 * t70 + t51 * t66;
t18 = -pkin(3) * t36 - pkin(10) * t37 + t26;
t65 = sin(qJ(4));
t9 = t82 * t14 + t65 * t18;
t15 = t70 * t23 - t66 * t27;
t35 = qJD(4) - t36;
t6 = t35 * qJ(5) + t9;
t78 = pkin(3) * t42 + t15;
t8 = -t65 * t14 + t82 * t18;
t75 = qJD(5) - t8;
t31 = t82 * t37 + t65 * t42;
t74 = qJ(5) * t31 + t78;
t73 = V_base(3) ^ 2;
t69 = cos(qJ(6));
t64 = sin(qJ(6));
t34 = qJD(6) - t35;
t30 = t37 * t65 - t82 * t42;
t20 = t30 * t64 + t31 * t69;
t19 = t30 * t69 - t31 * t64;
t10 = pkin(4) * t30 - t74;
t7 = t84 * t30 + t74;
t5 = -t35 * pkin(4) + t75;
t4 = pkin(11) * t30 + t6;
t3 = -t31 * pkin(11) + t84 * t35 + t75;
t2 = t3 * t64 + t4 * t69;
t1 = t3 * t69 - t4 * t64;
t11 = m(2) * (t52 ^ 2 + t53 ^ 2 + t73) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t73) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t42 / 0.2e1) * t42 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(3) * (t28 ^ 2 + t29 ^ 2 + t32 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t26 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(2,1) + t53 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,6) * t61 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t32 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,5) * t51 + Ifges(3,1) * t45 / 0.2e1) * t45 + (-t32 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t45 + Ifges(3,6) * t51 + Ifges(3,2) * t44 / 0.2e1) * t44 + (t52 * mrSges(2,1) - t53 * mrSges(2,2) + Ifges(2,3) * t61 / 0.2e1) * t61 + (-t26 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,6) * t42 + Ifges(4,2) * t36 / 0.2e1) * t36 + (t8 * mrSges(5,1) - t5 * mrSges(6,1) - t9 * mrSges(5,2) + t6 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t35) * t35 + (t26 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,5) * t42 + Ifges(4,1) * t37 / 0.2e1) * t37 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t34 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t20 / 0.2e1) * t20 + m(5) * (t78 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (-t78 * mrSges(5,1) + t10 * mrSges(6,1) - t6 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t30 + (-Ifges(5,6) + Ifges(6,6)) * t35 + (-Ifges(5,4) + Ifges(6,5)) * t31) * t30 + (-t78 * mrSges(5,2) + t5 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t31 + (Ifges(6,4) + Ifges(5,5)) * t35) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t28 * mrSges(3,1) - t29 * mrSges(3,2) + Ifges(3,3) * t51 / 0.2e1) * t51 + (V_base(3) * mrSges(2,2) - t52 * mrSges(2,3) + Ifges(2,5) * t61 + Ifges(2,1) * t55 / 0.2e1) * t55;
T  = t11;

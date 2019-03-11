% Calculate kinetic energy for
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:02
% EndTime: 2019-03-08 23:31:04
% DurationCPUTime: 1.36s
% Computational Cost: add. (3781->158), mult. (6449->222), div. (0->0), fcn. (5406->12), ass. (0->62)
t59 = V_base(5) * qJ(1) + V_base(1);
t60 = -V_base(4) * qJ(1) + V_base(2);
t62 = sin(pkin(11));
t64 = cos(pkin(11));
t52 = -t59 * t62 + t64 * t60;
t65 = cos(pkin(6));
t55 = t62 * V_base(5) + t64 * V_base(4);
t82 = pkin(7) * t55;
t46 = V_base(6) * pkin(1) - t65 * t82 + t52;
t54 = -t62 * V_base(4) + t64 * V_base(5);
t61 = V_base(3) + qJD(1);
t63 = sin(pkin(6));
t49 = -pkin(1) * t54 - t63 * t82 + t61;
t84 = t46 * t65 + t49 * t63;
t53 = t64 * t59 + t62 * t60;
t74 = t54 * t65 + t63 * V_base(6);
t42 = pkin(7) * t74 + t53;
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t28 = -t69 * t42 + t72 * t84;
t44 = -t69 * t55 + t72 * t74;
t83 = -pkin(4) - pkin(5);
t81 = cos(qJ(4));
t32 = -t46 * t63 + t65 * t49;
t45 = t55 * t72 + t69 * t74;
t23 = -pkin(2) * t44 - pkin(8) * t45 + t32;
t29 = t72 * t42 + t69 * t84;
t51 = -t54 * t63 + t65 * V_base(6) + qJD(2);
t27 = pkin(8) * t51 + t29;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t16 = t68 * t23 + t71 * t27;
t43 = qJD(3) - t44;
t14 = pkin(9) * t43 + t16;
t26 = -t51 * pkin(2) - t28;
t36 = -t45 * t68 + t71 * t51;
t37 = t45 * t71 + t51 * t68;
t18 = -t36 * pkin(3) - t37 * pkin(9) + t26;
t67 = sin(qJ(4));
t9 = t81 * t14 + t67 * t18;
t15 = t71 * t23 - t68 * t27;
t35 = qJD(4) - t36;
t6 = t35 * qJ(5) + t9;
t77 = pkin(3) * t43 + t15;
t8 = -t67 * t14 + t18 * t81;
t75 = qJD(5) - t8;
t31 = t37 * t81 + t67 * t43;
t73 = qJ(5) * t31 + t77;
t70 = cos(qJ(6));
t66 = sin(qJ(6));
t34 = qJD(6) - t35;
t30 = t37 * t67 - t43 * t81;
t20 = t30 * t66 + t31 * t70;
t19 = t30 * t70 - t31 * t66;
t10 = pkin(4) * t30 - t73;
t7 = t30 * t83 + t73;
t5 = -t35 * pkin(4) + t75;
t4 = pkin(10) * t30 + t6;
t3 = -t31 * pkin(10) + t35 * t83 + t75;
t2 = t3 * t66 + t4 * t70;
t1 = t3 * t70 - t4 * t66;
t11 = m(2) * (t52 ^ 2 + t53 ^ 2 + t61 ^ 2) / 0.2e1 + m(3) * (t28 ^ 2 + t29 ^ 2 + t32 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t61 * mrSges(2,2) - t52 * mrSges(2,3) + Ifges(2,1) * t55 / 0.2e1) * t55 + (t28 * mrSges(3,1) - t29 * mrSges(3,2) + Ifges(3,3) * t51 / 0.2e1) * t51 + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t43 / 0.2e1) * t43 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t61 * mrSges(2,1) + t53 * mrSges(2,3) + Ifges(2,4) * t55 + Ifges(2,2) * t54 / 0.2e1) * t54 + (t32 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,5) * t51 + Ifges(3,1) * t45 / 0.2e1) * t45 + (t26 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,5) * t43 + Ifges(4,1) * t37 / 0.2e1) * t37 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t20 / 0.2e1) * t20 + (-t32 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t45 + Ifges(3,6) * t51 + Ifges(3,2) * t44 / 0.2e1) * t44 + (-t26 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t37 + Ifges(4,6) * t43 + Ifges(4,2) * t36 / 0.2e1) * t36 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t34 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t8 * mrSges(5,1) - t5 * mrSges(6,1) - t9 * mrSges(5,2) + t6 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t35) * t35 + (-t77 * mrSges(5,2) + t5 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t31 + (Ifges(6,4) + Ifges(5,5)) * t35) * t31 + (V_base(2) * mrSges(1,1) + t52 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t53 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t55 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t54 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (-t77 * mrSges(5,1) + t10 * mrSges(6,1) - t6 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t30 + (-Ifges(5,6) + Ifges(6,6)) * t35 + (-Ifges(5,4) + Ifges(6,5)) * t31) * t30;
T  = t11;

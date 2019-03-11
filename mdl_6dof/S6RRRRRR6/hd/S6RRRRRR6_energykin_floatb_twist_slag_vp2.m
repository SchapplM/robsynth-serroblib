% Calculate kinetic energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:35
% EndTime: 2019-03-10 04:11:36
% DurationCPUTime: 1.65s
% Computational Cost: add. (6537->161), mult. (9743->240), div. (0->0), fcn. (8340->14), ass. (0->67)
t65 = pkin(7) * V_base(5) + V_base(1);
t66 = -pkin(7) * V_base(4) + V_base(2);
t75 = sin(qJ(1));
t81 = cos(qJ(1));
t57 = -t65 * t75 + t66 * t81;
t67 = V_base(6) + qJD(1);
t69 = cos(pkin(6));
t61 = t75 * V_base(5) + t81 * V_base(4);
t87 = pkin(8) * t61;
t52 = pkin(1) * t67 - t69 * t87 + t57;
t60 = -t75 * V_base(4) + t81 * V_base(5);
t68 = sin(pkin(6));
t55 = -pkin(1) * t60 - t68 * t87 + V_base(3);
t88 = t52 * t69 + t55 * t68;
t58 = t65 * t81 + t66 * t75;
t83 = t60 * t69 + t67 * t68;
t49 = pkin(8) * t83 + t58;
t74 = sin(qJ(2));
t80 = cos(qJ(2));
t38 = -t49 * t74 + t80 * t88;
t50 = -t61 * t74 + t80 * t83;
t40 = -t52 * t68 + t55 * t69;
t51 = t61 * t80 + t74 * t83;
t30 = -pkin(2) * t50 - pkin(9) * t51 + t40;
t39 = t80 * t49 + t74 * t88;
t56 = -t60 * t68 + t67 * t69 + qJD(2);
t34 = pkin(9) * t56 + t39;
t73 = sin(qJ(3));
t79 = cos(qJ(3));
t23 = t30 * t79 - t34 * t73;
t43 = t51 * t79 + t56 * t73;
t48 = qJD(3) - t50;
t17 = pkin(3) * t48 - pkin(10) * t43 + t23;
t24 = t30 * t73 + t34 * t79;
t42 = -t51 * t73 + t56 * t79;
t22 = pkin(10) * t42 + t24;
t72 = sin(qJ(4));
t78 = cos(qJ(4));
t12 = t17 * t72 + t22 * t78;
t47 = qJD(4) + t48;
t10 = pkin(11) * t47 + t12;
t33 = -pkin(2) * t56 - t38;
t25 = -pkin(3) * t42 + t33;
t36 = t42 * t78 - t43 * t72;
t37 = t42 * t72 + t43 * t78;
t15 = -pkin(4) * t36 - pkin(11) * t37 + t25;
t71 = sin(qJ(5));
t77 = cos(qJ(5));
t6 = t10 * t77 + t15 * t71;
t5 = -t10 * t71 + t15 * t77;
t11 = t17 * t78 - t22 * t72;
t35 = qJD(5) - t36;
t9 = -pkin(4) * t47 - t11;
t82 = V_base(3) ^ 2;
t76 = cos(qJ(6));
t70 = sin(qJ(6));
t32 = qJD(6) + t35;
t27 = t37 * t77 + t47 * t71;
t26 = -t37 * t71 + t47 * t77;
t19 = t26 * t70 + t27 * t76;
t18 = t26 * t76 - t27 * t70;
t7 = -pkin(5) * t26 + t9;
t4 = pkin(12) * t26 + t6;
t3 = pkin(5) * t35 - pkin(12) * t27 + t5;
t2 = t3 * t70 + t4 * t76;
t1 = t3 * t76 - t4 * t70;
t8 = (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t35 / 0.2e1) * t35 + (-t33 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t48 + Ifges(4,2) * t42 / 0.2e1) * t42 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t32 / 0.2e1) * t32 + (t57 * mrSges(2,1) - t58 * mrSges(2,2) + Ifges(2,3) * t67 / 0.2e1) * t67 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t32 + Ifges(7,1) * t19 / 0.2e1) * t19 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + (-t25 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,6) * t47 + Ifges(5,2) * t36 / 0.2e1) * t36 + (-t40 * mrSges(3,1) + t39 * mrSges(3,3) + Ifges(3,4) * t51 + Ifges(3,6) * t56 + Ifges(3,2) * t50 / 0.2e1) * t50 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t35 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t35 + Ifges(6,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t25 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t47 + Ifges(5,1) * t37 / 0.2e1) * t37 + (t33 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t43 / 0.2e1) * t43 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t32 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t23 * mrSges(4,1) - t24 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (t38 * mrSges(3,1) - t39 * mrSges(3,2) + Ifges(3,3) * t56 / 0.2e1) * t56 + (t40 * mrSges(3,2) - t38 * mrSges(3,3) + Ifges(3,5) * t56 + Ifges(3,1) * t51 / 0.2e1) * t51 + (V_base(3) * mrSges(2,2) - t57 * mrSges(2,3) + Ifges(2,5) * t67 + Ifges(2,1) * t61 / 0.2e1) * t61 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(2,1) + t58 * mrSges(2,3) + Ifges(2,4) * t61 + Ifges(2,6) * t67 + Ifges(2,2) * t60 / 0.2e1) * t60 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t33 ^ 2) / 0.2e1 + m(3) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t82) / 0.2e1 + m(2) * (t57 ^ 2 + t58 ^ 2 + t82) / 0.2e1;
T  = t8;

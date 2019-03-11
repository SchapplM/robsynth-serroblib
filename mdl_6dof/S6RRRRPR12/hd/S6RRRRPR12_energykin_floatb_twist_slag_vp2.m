% Calculate kinetic energy for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:13
% EndTime: 2019-03-09 23:29:15
% DurationCPUTime: 1.73s
% Computational Cost: add. (10447->166), mult. (16419->248), div. (0->0), fcn. (14354->16), ass. (0->72)
t65 = pkin(8) * V_base(5) + V_base(1);
t66 = -pkin(8) * V_base(4) + V_base(2);
t78 = sin(qJ(1));
t83 = cos(qJ(1));
t58 = -t65 * t78 + t83 * t66;
t67 = V_base(6) + qJD(1);
t73 = cos(pkin(6));
t61 = t78 * V_base(5) + t83 * V_base(4);
t94 = pkin(9) * t61;
t52 = pkin(1) * t67 - t73 * t94 + t58;
t60 = -t78 * V_base(4) + t83 * V_base(5);
t70 = sin(pkin(6));
t56 = -pkin(1) * t60 - t70 * t94 + V_base(3);
t96 = t52 * t73 + t56 * t70;
t59 = t83 * t65 + t78 * t66;
t85 = t60 * t73 + t67 * t70;
t49 = t85 * pkin(9) + t59;
t77 = sin(qJ(2));
t82 = cos(qJ(2));
t39 = -t49 * t77 + t96 * t82;
t57 = -t60 * t70 + t67 * t73 + qJD(2);
t72 = cos(pkin(7));
t51 = t61 * t82 + t85 * t77;
t93 = pkin(10) * t51;
t35 = pkin(2) * t57 - t72 * t93 + t39;
t44 = -t52 * t70 + t73 * t56;
t50 = -t61 * t77 + t85 * t82;
t69 = sin(pkin(7));
t38 = -pkin(2) * t50 - t69 * t93 + t44;
t95 = t35 * t72 + t38 * t69;
t40 = t82 * t49 + t96 * t77;
t86 = t50 * t72 + t57 * t69;
t34 = t86 * pkin(10) + t40;
t76 = sin(qJ(3));
t81 = cos(qJ(3));
t21 = -t76 * t34 + t95 * t81;
t42 = -t51 * t76 + t86 * t81;
t22 = t81 * t34 + t95 * t76;
t45 = -t50 * t69 + t57 * t72 + qJD(3);
t17 = pkin(11) * t45 + t22;
t28 = -t35 * t69 + t72 * t38;
t43 = t51 * t81 + t86 * t76;
t20 = -pkin(3) * t42 - pkin(11) * t43 + t28;
t75 = sin(qJ(4));
t80 = cos(qJ(4));
t13 = t80 * t17 + t75 * t20;
t32 = -t43 * t75 + t45 * t80;
t11 = qJ(5) * t32 + t13;
t68 = sin(pkin(13));
t71 = cos(pkin(13));
t12 = -t17 * t75 + t80 * t20;
t33 = t43 * t80 + t45 * t75;
t41 = qJD(4) - t42;
t9 = pkin(4) * t41 - qJ(5) * t33 + t12;
t6 = t71 * t11 + t68 * t9;
t5 = -t11 * t68 + t71 * t9;
t26 = t32 * t71 - t33 * t68;
t16 = -pkin(3) * t45 - t21;
t14 = -pkin(4) * t32 + qJD(5) + t16;
t84 = V_base(3) ^ 2;
t79 = cos(qJ(6));
t74 = sin(qJ(6));
t27 = t32 * t68 + t33 * t71;
t25 = qJD(6) - t26;
t24 = t27 * t79 + t41 * t74;
t23 = -t27 * t74 + t41 * t79;
t7 = -pkin(5) * t26 - pkin(12) * t27 + t14;
t4 = pkin(12) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t4 * t79 + t7 * t74;
t1 = -t4 * t74 + t7 * t79;
t8 = (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t24 / 0.2e1) * t24 + (-t16 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,2) * t32 / 0.2e1) * t32 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t24 + Ifges(7,6) * t25 + Ifges(7,2) * t23 / 0.2e1) * t23 + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t13 * mrSges(5,2) - t6 * mrSges(6,2) + Ifges(5,5) * t33 + Ifges(6,5) * t27 + Ifges(5,6) * t32 + Ifges(6,6) * t26 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t41) * t41 + (V_base(3) * mrSges(2,2) - t58 * mrSges(2,3) + Ifges(2,5) * t67 + Ifges(2,1) * t61 / 0.2e1) * t61 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t27 / 0.2e1) * t27 + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,2) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t44 * mrSges(3,2) - t39 * mrSges(3,3) + Ifges(3,5) * t57 + Ifges(3,1) * t51 / 0.2e1) * t51 + (t28 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t45 + Ifges(4,1) * t43 / 0.2e1) * t43 + (-t28 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t45 + Ifges(4,2) * t42 / 0.2e1) * t42 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (t58 * mrSges(2,1) - t59 * mrSges(2,2) + Ifges(2,3) * t67 / 0.2e1) * t67 + (-V_base(3) * mrSges(2,1) + t59 * mrSges(2,3) + Ifges(2,4) * t61 + Ifges(2,6) * t67 + Ifges(2,2) * t60 / 0.2e1) * t60 + (t39 * mrSges(3,1) - t40 * mrSges(3,2) + Ifges(3,3) * t57 / 0.2e1) * t57 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t45 / 0.2e1) * t45 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (-t44 * mrSges(3,1) + t40 * mrSges(3,3) + Ifges(3,4) * t51 + Ifges(3,6) * t57 + Ifges(3,2) * t50 / 0.2e1) * t50 + (t16 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t33 / 0.2e1) * t33 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t28 ^ 2) / 0.2e1 + m(3) * (t39 ^ 2 + t40 ^ 2 + t44 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t84) / 0.2e1 + m(2) * (t58 ^ 2 + t59 ^ 2 + t84) / 0.2e1;
T  = t8;

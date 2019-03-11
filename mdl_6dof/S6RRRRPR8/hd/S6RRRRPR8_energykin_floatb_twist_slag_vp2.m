% Calculate kinetic energy for
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:44
% EndTime: 2019-03-09 22:38:45
% DurationCPUTime: 1.29s
% Computational Cost: add. (2877->153), mult. (3645->213), div. (0->0), fcn. (2876->10), ass. (0->57)
t75 = -pkin(4) - pkin(5);
t74 = cos(qJ(4));
t64 = sin(qJ(1));
t68 = cos(qJ(1));
t51 = -t64 * V_base(4) + t68 * V_base(5);
t52 = t64 * V_base(5) + t68 * V_base(4);
t36 = -pkin(1) * t51 - pkin(7) * t52 + V_base(3);
t57 = V_base(5) * pkin(6) + V_base(1);
t58 = -V_base(4) * pkin(6) + V_base(2);
t48 = t68 * t57 + t64 * t58;
t59 = V_base(6) + qJD(1);
t43 = pkin(7) * t59 + t48;
t63 = sin(qJ(2));
t67 = cos(qJ(2));
t30 = t63 * t36 + t67 * t43;
t50 = qJD(2) - t51;
t27 = pkin(8) * t50 + t30;
t47 = -t64 * t57 + t58 * t68;
t42 = -pkin(1) * t59 - t47;
t45 = -t63 * t52 + t67 * t59;
t46 = t52 * t67 + t59 * t63;
t28 = -pkin(2) * t45 - pkin(8) * t46 + t42;
t62 = sin(qJ(3));
t66 = cos(qJ(3));
t18 = -t27 * t62 + t66 * t28;
t33 = t46 * t66 + t50 * t62;
t44 = qJD(3) - t45;
t12 = pkin(3) * t44 - pkin(9) * t33 + t18;
t19 = t66 * t27 + t62 * t28;
t32 = -t46 * t62 + t50 * t66;
t15 = pkin(9) * t32 + t19;
t61 = sin(qJ(4));
t8 = t61 * t12 + t74 * t15;
t29 = t67 * t36 - t63 * t43;
t41 = qJD(4) + t44;
t6 = t41 * qJ(5) + t8;
t73 = pkin(2) * t50 + t29;
t7 = t12 * t74 - t61 * t15;
t72 = qJD(5) - t7;
t71 = pkin(3) * t32 + t73;
t22 = t61 * t32 + t33 * t74;
t70 = qJ(5) * t22 + t71;
t69 = V_base(3) ^ 2;
t65 = cos(qJ(6));
t60 = sin(qJ(6));
t40 = qJD(6) - t41;
t21 = -t32 * t74 + t33 * t61;
t17 = t21 * t60 + t22 * t65;
t16 = t21 * t65 - t22 * t60;
t10 = pkin(4) * t21 - t70;
t9 = t21 * t75 + t70;
t5 = -t41 * pkin(4) + t72;
t4 = pkin(10) * t21 + t6;
t3 = -t22 * pkin(10) + t41 * t75 + t72;
t2 = t3 * t60 + t4 * t65;
t1 = t3 * t65 - t4 * t60;
t11 = (t7 * mrSges(5,1) - t5 * mrSges(6,1) - t8 * mrSges(5,2) + t6 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t41) * t41 + (-t42 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t46 + Ifges(3,6) * t50 + Ifges(3,2) * t45 / 0.2e1) * t45 + (-t9 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t40 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t50 / 0.2e1) * t50 + (t9 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t40 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(2,2) - t47 * mrSges(2,3) + Ifges(2,5) * t59 + Ifges(2,1) * t52 / 0.2e1) * t52 + (-V_base(3) * mrSges(2,1) + t48 * mrSges(2,3) + Ifges(2,4) * t52 + Ifges(2,6) * t59 + Ifges(2,2) * t51 / 0.2e1) * t51 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t40 / 0.2e1) * t40 + (t18 * mrSges(4,1) - t19 * mrSges(4,2) + Ifges(4,3) * t44 / 0.2e1) * t44 + (t47 * mrSges(2,1) - t48 * mrSges(2,2) + Ifges(2,3) * t59 / 0.2e1) * t59 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(2) * (t47 ^ 2 + t48 ^ 2 + t69) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t69) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t42 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (-t71 * mrSges(5,2) + t5 * mrSges(6,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t22 + (Ifges(6,4) + Ifges(5,5)) * t41) * t22 + (-t71 * mrSges(5,1) + t10 * mrSges(6,1) - t6 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t21 + (-Ifges(5,6) + Ifges(6,6)) * t41 + (-Ifges(5,4) + Ifges(6,5)) * t22) * t21 + m(5) * (t7 ^ 2 + t71 ^ 2 + t8 ^ 2) / 0.2e1 + (-t73 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,5) * t44 + Ifges(4,1) * t33 / 0.2e1) * t33 + (t73 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t33 + Ifges(4,6) * t44 + Ifges(4,2) * t32 / 0.2e1) * t32 + m(4) * (t18 ^ 2 + t19 ^ 2 + t73 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t42 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t50 + Ifges(3,1) * t46 / 0.2e1) * t46;
T  = t11;

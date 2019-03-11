% Calculate kinetic energy for
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:31
% EndTime: 2019-03-09 19:15:32
% DurationCPUTime: 1.10s
% Computational Cost: add. (2633->153), mult. (3311->213), div. (0->0), fcn. (2580->10), ass. (0->55)
t66 = sin(qJ(1));
t70 = cos(qJ(1));
t54 = t66 * V_base(5) + t70 * V_base(4);
t61 = V_base(6) + qJD(1);
t65 = sin(qJ(2));
t69 = cos(qJ(2));
t47 = t54 * t69 + t61 * t65;
t53 = -t66 * V_base(4) + t70 * V_base(5);
t52 = qJD(2) - t53;
t64 = sin(qJ(3));
t73 = cos(qJ(3));
t34 = t47 * t73 + t64 * t52;
t46 = -t54 * t65 + t69 * t61;
t45 = qJD(3) - t46;
t37 = -pkin(1) * t53 - pkin(7) * t54 + V_base(3);
t59 = V_base(5) * pkin(6) + V_base(1);
t60 = -V_base(4) * pkin(6) + V_base(2);
t49 = t70 * t59 + t66 * t60;
t44 = pkin(7) * t61 + t49;
t30 = t65 * t37 + t69 * t44;
t27 = pkin(8) * t52 + t30;
t48 = -t66 * t59 + t60 * t70;
t43 = -pkin(1) * t61 - t48;
t28 = -pkin(2) * t46 - pkin(8) * t47 + t43;
t19 = -t64 * t27 + t28 * t73;
t72 = qJD(4) - t19;
t10 = -t34 * pkin(9) + (-pkin(3) - pkin(4)) * t45 + t72;
t20 = t73 * t27 + t64 * t28;
t17 = t45 * qJ(4) + t20;
t33 = t47 * t64 - t52 * t73;
t12 = pkin(9) * t33 + t17;
t63 = sin(qJ(5));
t68 = cos(qJ(5));
t6 = t63 * t10 + t68 * t12;
t29 = t69 * t37 - t65 * t44;
t26 = -t52 * pkin(2) - t29;
t5 = t68 * t10 - t12 * t63;
t18 = t33 * pkin(3) - t34 * qJ(4) + t26;
t42 = qJD(5) - t45;
t13 = -pkin(4) * t33 - t18;
t71 = V_base(3) ^ 2;
t67 = cos(qJ(6));
t62 = sin(qJ(6));
t40 = qJD(6) + t42;
t22 = t33 * t63 + t34 * t68;
t21 = t33 * t68 - t34 * t63;
t16 = -t45 * pkin(3) + t72;
t15 = t21 * t62 + t22 * t67;
t14 = t21 * t67 - t22 * t62;
t7 = -pkin(5) * t21 + t13;
t4 = pkin(10) * t21 + t6;
t3 = pkin(5) * t42 - pkin(10) * t22 + t5;
t2 = t3 * t62 + t4 * t67;
t1 = t3 * t67 - t4 * t62;
t8 = (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t19 * mrSges(4,1) - t16 * mrSges(5,1) - t20 * mrSges(4,2) + t17 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t45) * t45 + (t13 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t42 + Ifges(6,1) * t22 / 0.2e1) * t22 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t15 + Ifges(7,6) * t40 + Ifges(7,2) * t14 / 0.2e1) * t14 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t40 / 0.2e1) * t40 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t43 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t47 + Ifges(3,6) * t52 + Ifges(3,2) * t46 / 0.2e1) * t46 + (t26 * mrSges(4,2) + t16 * mrSges(5,2) - t19 * mrSges(4,3) - t18 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t34 + (Ifges(5,4) + Ifges(4,5)) * t45) * t34 + (t43 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t52 + Ifges(3,1) * t47 / 0.2e1) * t47 + (V_base(3) * mrSges(2,2) - t48 * mrSges(2,3) + Ifges(2,5) * t61 + Ifges(2,1) * t54 / 0.2e1) * t54 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t40 + Ifges(7,1) * t15 / 0.2e1) * t15 + m(2) * (t48 ^ 2 + t49 ^ 2 + t71) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t71) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t43 ^ 2) / 0.2e1 + m(4) * (t19 ^ 2 + t20 ^ 2 + t26 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t13 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + (t48 * mrSges(2,1) - t49 * mrSges(2,2) + Ifges(2,3) * t61 / 0.2e1) * t61 + (-V_base(3) * mrSges(2,1) + t49 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,6) * t61 + Ifges(2,2) * t53 / 0.2e1) * t53 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t52 / 0.2e1) * t52 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t42 / 0.2e1) * t42 + (-t13 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t22 + Ifges(6,6) * t42 + Ifges(6,2) * t21 / 0.2e1) * t21 + (t26 * mrSges(4,1) + t18 * mrSges(5,1) - t17 * mrSges(5,2) - t20 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t33 + (-Ifges(4,6) + Ifges(5,6)) * t45 + (-Ifges(4,4) + Ifges(5,5)) * t34) * t33;
T  = t8;

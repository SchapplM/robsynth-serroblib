% Calculate kinetic energy for
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR11_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:38:45
% EndTime: 2019-03-09 09:38:47
% DurationCPUTime: 1.14s
% Computational Cost: add. (3949->160), mult. (5927->225), div. (0->0), fcn. (4930->12), ass. (0->63)
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t51 = t66 * V_base(5) + t69 * V_base(4);
t78 = pkin(8) * t51;
t65 = sin(qJ(2));
t50 = -t66 * V_base(4) + t69 * V_base(5);
t58 = V_base(6) + qJD(1);
t60 = sin(pkin(6));
t62 = cos(pkin(6));
t73 = t50 * t62 + t58 * t60;
t77 = cos(qJ(2));
t40 = t77 * t51 + t73 * t65;
t46 = -t50 * t60 + t58 * t62 + qJD(2);
t56 = V_base(5) * pkin(7) + V_base(1);
t57 = -V_base(4) * pkin(7) + V_base(2);
t48 = t69 * t56 + t66 * t57;
t38 = t73 * pkin(8) + t48;
t47 = -t56 * t66 + t69 * t57;
t41 = pkin(1) * t58 - t62 * t78 + t47;
t44 = -pkin(1) * t50 - t60 * t78 + V_base(3);
t74 = t62 * t77;
t75 = t60 * t77;
t29 = -t65 * t38 + t41 * t74 + t44 * t75;
t71 = qJD(3) - t29;
t76 = pkin(2) + qJ(4);
t17 = t40 * pkin(3) - t76 * t46 + t71;
t39 = -t50 * t74 + t51 * t65 - t58 * t75;
t31 = -t41 * t60 + t62 * t44;
t72 = -qJ(3) * t40 + t31;
t19 = t76 * t39 + t72;
t59 = sin(pkin(11));
t61 = cos(pkin(11));
t13 = t59 * t17 + t61 * t19;
t32 = t39 * t61 - t46 * t59;
t11 = pkin(9) * t32 + t13;
t64 = sin(qJ(5));
t68 = cos(qJ(5));
t12 = t61 * t17 - t19 * t59;
t33 = t39 * t59 + t46 * t61;
t9 = pkin(4) * t40 - pkin(9) * t33 + t12;
t6 = t68 * t11 + t64 * t9;
t30 = t77 * t38 + (t41 * t62 + t44 * t60) * t65;
t26 = -t46 * qJ(3) - t30;
t5 = -t11 * t64 + t68 * t9;
t27 = t32 * t68 - t33 * t64;
t22 = -pkin(3) * t39 + qJD(4) - t26;
t14 = -pkin(4) * t32 + t22;
t70 = V_base(3) ^ 2;
t67 = cos(qJ(6));
t63 = sin(qJ(6));
t37 = qJD(5) + t40;
t28 = t32 * t64 + t33 * t68;
t25 = qJD(6) - t27;
t24 = -t46 * pkin(2) + t71;
t23 = pkin(2) * t39 + t72;
t21 = t28 * t67 + t37 * t63;
t20 = -t28 * t63 + t37 * t67;
t7 = -pkin(5) * t27 - pkin(10) * t28 + t14;
t4 = pkin(10) * t37 + t6;
t3 = -pkin(5) * t37 - t5;
t2 = t4 * t67 + t63 * t7;
t1 = -t4 * t63 + t67 * t7;
t8 = m(2) * (t47 ^ 2 + t48 ^ 2 + t70) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t70) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t22 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t26 ^ 2) / 0.2e1 + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t47 * mrSges(2,1) - t48 * mrSges(2,2) + Ifges(2,3) * t58 / 0.2e1) * t58 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t37 / 0.2e1) * t37 + (t22 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t33 / 0.2e1) * t33 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t47 * mrSges(2,3) + Ifges(2,5) * t58 + Ifges(2,1) * t51 / 0.2e1) * t51 + (-t22 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,2) * t32 / 0.2e1) * t32 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t37 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t21 / 0.2e1) * t21 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t48 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,6) * t58 + Ifges(2,2) * t50 / 0.2e1) * t50 + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t37 + Ifges(6,2) * t27 / 0.2e1) * t27 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t25 + Ifges(7,2) * t20 / 0.2e1) * t20 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + t24 * mrSges(4,2) - t26 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t46) * t46 + (t31 * mrSges(3,1) + t26 * mrSges(4,1) - t23 * mrSges(4,2) - t30 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t39 + (Ifges(4,5) - Ifges(3,6)) * t46) * t39 + (t24 * mrSges(4,1) + t12 * mrSges(5,1) + t31 * mrSges(3,2) - t13 * mrSges(5,2) - t29 * mrSges(3,3) - t23 * mrSges(4,3) + Ifges(5,5) * t33 + Ifges(5,6) * t32 + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t40 + (-Ifges(4,4) + Ifges(3,5)) * t46 + (-Ifges(3,4) - Ifges(4,6)) * t39) * t40;
T  = t8;

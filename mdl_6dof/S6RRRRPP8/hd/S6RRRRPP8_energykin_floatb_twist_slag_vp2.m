% Calculate kinetic energy for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:12
% EndTime: 2019-03-09 21:30:13
% DurationCPUTime: 1.16s
% Computational Cost: add. (3539->154), mult. (5277->206), div. (0->0), fcn. (4376->10), ass. (0->56)
t54 = V_base(5) * pkin(7) + V_base(1);
t55 = -V_base(4) * pkin(7) + V_base(2);
t62 = sin(qJ(1));
t64 = cos(qJ(1));
t47 = -t54 * t62 + t64 * t55;
t56 = V_base(6) + qJD(1);
t58 = cos(pkin(6));
t50 = t62 * V_base(5) + t64 * V_base(4);
t76 = pkin(8) * t50;
t41 = pkin(1) * t56 - t58 * t76 + t47;
t49 = -t62 * V_base(4) + t64 * V_base(5);
t57 = sin(pkin(6));
t44 = -pkin(1) * t49 - t57 * t76 + V_base(3);
t78 = t41 * t58 + t44 * t57;
t48 = t64 * t54 + t62 * t55;
t68 = t49 * t58 + t56 * t57;
t38 = pkin(8) * t68 + t48;
t61 = sin(qJ(2));
t63 = cos(qJ(2));
t24 = -t61 * t38 + t63 * t78;
t39 = -t50 * t61 + t63 * t68;
t77 = -pkin(4) - pkin(5);
t75 = cos(qJ(3));
t74 = cos(qJ(4));
t28 = -t41 * t57 + t58 * t44;
t40 = t50 * t63 + t61 * t68;
t19 = -pkin(2) * t39 - pkin(9) * t40 + t28;
t25 = t63 * t38 + t61 * t78;
t46 = -t49 * t57 + t56 * t58 + qJD(2);
t23 = pkin(9) * t46 + t25;
t60 = sin(qJ(3));
t14 = t60 * t19 + t75 * t23;
t37 = qJD(3) - t39;
t12 = pkin(10) * t37 + t14;
t22 = -pkin(2) * t46 - t24;
t31 = -t40 * t60 + t75 * t46;
t32 = t40 * t75 + t60 * t46;
t16 = -pkin(3) * t31 - pkin(10) * t32 + t22;
t59 = sin(qJ(4));
t7 = t74 * t12 + t59 * t16;
t13 = t75 * t19 - t60 * t23;
t30 = qJD(4) - t31;
t4 = t30 * qJ(5) + t7;
t70 = pkin(3) * t37 + t13;
t6 = -t59 * t12 + t16 * t74;
t67 = qJD(5) - t6;
t27 = t32 * t74 + t59 * t37;
t66 = qJ(5) * t27 + t70;
t65 = V_base(3) ^ 2;
t26 = t32 * t59 - t37 * t74;
t8 = pkin(4) * t26 - t66;
t5 = t26 * t77 + qJD(6) + t66;
t3 = -t30 * pkin(4) + t67;
t2 = qJ(6) * t26 + t4;
t1 = -t27 * qJ(6) + t77 * t30 + t67;
t9 = m(2) * (t47 ^ 2 + t48 ^ 2 + t65) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t65) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t22 ^ 2) / 0.2e1 + m(5) * (t6 ^ 2 + t7 ^ 2 + t70 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t47 * mrSges(2,1) - t48 * mrSges(2,2) + Ifges(2,3) * t56 / 0.2e1) * t56 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t37 / 0.2e1) * t37 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t47 * mrSges(2,3) + Ifges(2,5) * t56 + Ifges(2,1) * t50 / 0.2e1) * t50 + (t28 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t40 / 0.2e1) * t40 + (t22 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t37 + Ifges(4,1) * t32 / 0.2e1) * t32 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t48 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,6) * t56 + Ifges(2,2) * t49 / 0.2e1) * t49 + (-t28 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t40 + Ifges(3,6) * t46 + Ifges(3,2) * t39 / 0.2e1) * t39 + (-t22 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t37 + Ifges(4,2) * t31 / 0.2e1) * t31 + (t6 * mrSges(5,1) - t3 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t30) * t30 + (-t70 * mrSges(5,2) + t3 * mrSges(6,2) + t5 * mrSges(7,2) - t6 * mrSges(5,3) - t8 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t27 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t30) * t27 + (-t70 * mrSges(5,1) + t8 * mrSges(6,1) - t5 * mrSges(7,1) - t4 * mrSges(6,2) - t7 * mrSges(5,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t26 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t30 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t27) * t26;
T  = t9;

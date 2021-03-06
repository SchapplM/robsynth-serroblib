% Calculate kinetic energy for
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:38
% EndTime: 2019-03-08 22:48:39
% DurationCPUTime: 1.11s
% Computational Cost: add. (3095->154), mult. (5277->205), div. (0->0), fcn. (4376->10), ass. (0->55)
t54 = V_base(5) * qJ(1) + V_base(1);
t55 = -V_base(4) * qJ(1) + V_base(2);
t57 = sin(pkin(10));
t59 = cos(pkin(10));
t47 = -t54 * t57 + t59 * t55;
t60 = cos(pkin(6));
t50 = t57 * V_base(5) + t59 * V_base(4);
t75 = pkin(7) * t50;
t41 = V_base(6) * pkin(1) - t60 * t75 + t47;
t49 = -t57 * V_base(4) + t59 * V_base(5);
t56 = V_base(3) + qJD(1);
t58 = sin(pkin(6));
t44 = -pkin(1) * t49 - t58 * t75 + t56;
t77 = t41 * t60 + t44 * t58;
t48 = t59 * t54 + t57 * t55;
t66 = t49 * t60 + t58 * V_base(6);
t37 = t66 * pkin(7) + t48;
t63 = sin(qJ(2));
t64 = cos(qJ(2));
t24 = -t63 * t37 + t64 * t77;
t39 = -t63 * t50 + t66 * t64;
t76 = -pkin(4) - pkin(5);
t74 = cos(qJ(3));
t73 = cos(qJ(4));
t28 = -t41 * t58 + t60 * t44;
t40 = t50 * t64 + t66 * t63;
t19 = -pkin(2) * t39 - pkin(8) * t40 + t28;
t25 = t64 * t37 + t63 * t77;
t46 = -t49 * t58 + t60 * V_base(6) + qJD(2);
t23 = pkin(8) * t46 + t25;
t62 = sin(qJ(3));
t14 = t62 * t19 + t74 * t23;
t38 = qJD(3) - t39;
t12 = pkin(9) * t38 + t14;
t22 = -t46 * pkin(2) - t24;
t31 = -t40 * t62 + t74 * t46;
t32 = t74 * t40 + t62 * t46;
t16 = -t31 * pkin(3) - t32 * pkin(9) + t22;
t61 = sin(qJ(4));
t7 = t73 * t12 + t61 * t16;
t13 = t74 * t19 - t62 * t23;
t30 = qJD(4) - t31;
t4 = t30 * qJ(5) + t7;
t69 = pkin(3) * t38 + t13;
t6 = -t61 * t12 + t73 * t16;
t67 = qJD(5) - t6;
t27 = t73 * t32 + t61 * t38;
t65 = qJ(5) * t27 + t69;
t26 = t32 * t61 - t73 * t38;
t8 = pkin(4) * t26 - t65;
t5 = t76 * t26 + qJD(6) + t65;
t3 = -t30 * pkin(4) + t67;
t2 = qJ(6) * t26 + t4;
t1 = -t27 * qJ(6) + t76 * t30 + t67;
t9 = m(2) * (t47 ^ 2 + t48 ^ 2 + t56 ^ 2) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t22 ^ 2) / 0.2e1 + m(5) * (t6 ^ 2 + t69 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t56 * mrSges(2,2) - t47 * mrSges(2,3) + Ifges(2,1) * t50 / 0.2e1) * t50 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t56 * mrSges(2,1) + t48 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,2) * t49 / 0.2e1) * t49 + (t28 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t40 / 0.2e1) * t40 + (t22 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t38 + Ifges(4,1) * t32 / 0.2e1) * t32 + (-t28 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t40 + Ifges(3,6) * t46 + Ifges(3,2) * t39 / 0.2e1) * t39 + (-t22 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t38 + Ifges(4,2) * t31 / 0.2e1) * t31 + (V_base(2) * mrSges(1,1) + t47 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t48 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t50 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t49 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t6 * mrSges(5,1) - t3 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t30) * t30 + (-t69 * mrSges(5,2) + t3 * mrSges(6,2) + t5 * mrSges(7,2) - t6 * mrSges(5,3) - t8 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t27 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t30) * t27 + (-t69 * mrSges(5,1) + t8 * mrSges(6,1) - t5 * mrSges(7,1) - t4 * mrSges(6,2) - t7 * mrSges(5,3) + t2 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t26 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t30 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t27) * t26;
T  = t9;

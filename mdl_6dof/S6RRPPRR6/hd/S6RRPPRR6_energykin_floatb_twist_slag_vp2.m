% Calculate kinetic energy for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:12:56
% EndTime: 2019-03-09 09:12:57
% DurationCPUTime: 1.01s
% Computational Cost: add. (2483->153), mult. (3205->211), div. (0->0), fcn. (2476->10), ass. (0->54)
t65 = sin(qJ(1));
t71 = cos(qJ(1));
t49 = t65 * V_base(5) + t71 * V_base(4);
t59 = V_base(6) + qJD(1);
t64 = sin(qJ(2));
t70 = cos(qJ(2));
t42 = t70 * t49 + t64 * t59;
t48 = -t65 * V_base(4) + t71 * V_base(5);
t47 = qJD(2) - t48;
t34 = -pkin(1) * t48 - pkin(7) * t49 + V_base(3);
t55 = V_base(5) * pkin(6) + V_base(1);
t56 = -V_base(4) * pkin(6) + V_base(2);
t44 = t71 * t55 + t65 * t56;
t38 = pkin(7) * t59 + t44;
t29 = t70 * t34 - t64 * t38;
t69 = qJD(3) - t29;
t19 = -t42 * qJ(4) + (-pkin(2) - pkin(3)) * t47 + t69;
t30 = t64 * t34 + t70 * t38;
t27 = t47 * qJ(3) + t30;
t41 = t49 * t64 - t70 * t59;
t24 = qJ(4) * t41 + t27;
t60 = sin(pkin(10));
t61 = cos(pkin(10));
t13 = t60 * t19 + t61 * t24;
t31 = t41 * t61 - t42 * t60;
t11 = pkin(8) * t31 + t13;
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t12 = t61 * t19 - t24 * t60;
t32 = t41 * t60 + t42 * t61;
t9 = -pkin(4) * t47 - pkin(8) * t32 + t12;
t6 = t67 * t11 + t63 * t9;
t43 = -t65 * t55 + t71 * t56;
t37 = -t59 * pkin(1) - t43;
t28 = t41 * pkin(2) - t42 * qJ(3) + t37;
t5 = -t11 * t63 + t67 * t9;
t21 = t31 * t67 - t32 * t63;
t25 = -pkin(3) * t41 + qJD(4) - t28;
t14 = -pkin(4) * t31 + t25;
t68 = V_base(3) ^ 2;
t66 = cos(qJ(6));
t62 = sin(qJ(6));
t46 = qJD(5) - t47;
t26 = -t47 * pkin(2) + t69;
t22 = t31 * t63 + t32 * t67;
t20 = qJD(6) - t21;
t16 = t22 * t66 + t46 * t62;
t15 = -t22 * t62 + t46 * t66;
t7 = -pkin(5) * t21 - pkin(9) * t22 + t14;
t4 = pkin(9) * t46 + t6;
t3 = -pkin(5) * t46 - t5;
t2 = t4 * t66 + t62 * t7;
t1 = -t4 * t62 + t66 * t7;
t8 = m(2) * (t43 ^ 2 + t44 ^ 2 + t68) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t68) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t37 ^ 2) / 0.2e1 + m(4) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t25 ^ 2) / 0.2e1 + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + Ifges(2,3) * t59 / 0.2e1) * t59 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t46 / 0.2e1) * t46 + (t25 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t32 / 0.2e1) * t32 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,5) * t59 + Ifges(2,1) * t49 / 0.2e1) * t49 + (-t25 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t32 + Ifges(5,2) * t31 / 0.2e1) * t31 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t46 + Ifges(6,1) * t22 / 0.2e1) * t22 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t20 + Ifges(7,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t49 + Ifges(2,6) * t59 + Ifges(2,2) * t48 / 0.2e1) * t48 + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t22 + Ifges(6,6) * t46 + Ifges(6,2) * t21 / 0.2e1) * t21 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t20 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t37 * mrSges(3,2) + t26 * mrSges(4,2) - t29 * mrSges(3,3) - t28 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t42) * t42 + (t37 * mrSges(3,1) + t28 * mrSges(4,1) - t27 * mrSges(4,2) - t30 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t41 + (-Ifges(3,4) + Ifges(4,5)) * t42) * t41 + (t29 * mrSges(3,1) - t26 * mrSges(4,1) - t12 * mrSges(5,1) - t30 * mrSges(3,2) + t13 * mrSges(5,2) + t27 * mrSges(4,3) - Ifges(5,5) * t32 - Ifges(5,6) * t31 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t47 + (Ifges(4,4) + Ifges(3,5)) * t42 + (-Ifges(3,6) + Ifges(4,6)) * t41) * t47;
T  = t8;

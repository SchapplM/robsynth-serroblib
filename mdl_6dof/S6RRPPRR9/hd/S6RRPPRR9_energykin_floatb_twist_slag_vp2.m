% Calculate kinetic energy for
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:13
% EndTime: 2019-03-09 09:28:14
% DurationCPUTime: 0.94s
% Computational Cost: add. (2477->158), mult. (3717->210), div. (0->0), fcn. (2988->10), ass. (0->57)
t63 = sin(qJ(1));
t66 = cos(qJ(1));
t50 = t63 * V_base(5) + t66 * V_base(4);
t74 = pkin(8) * t50;
t49 = -t63 * V_base(4) + t66 * V_base(5);
t57 = V_base(6) + qJD(1);
t62 = sin(qJ(2));
t59 = cos(pkin(6));
t73 = cos(qJ(2));
t71 = t59 * t73;
t58 = sin(pkin(6));
t72 = t58 * t73;
t35 = -t49 * t71 + t50 * t62 - t57 * t72;
t45 = -t49 * t58 + t57 * t59 + qJD(2);
t55 = V_base(5) * pkin(7) + V_base(1);
t56 = -V_base(4) * pkin(7) + V_base(2);
t47 = t66 * t55 + t63 * t56;
t69 = t49 * t59 + t57 * t58;
t34 = pkin(8) * t69 + t47;
t46 = -t55 * t63 + t66 * t56;
t37 = pkin(1) * t57 - t59 * t74 + t46;
t41 = -pkin(1) * t49 - t58 * t74 + V_base(3);
t20 = t73 * t34 + (t37 * t59 + t41 * t58) * t62;
t18 = -t45 * qJ(3) - t20;
t70 = qJD(4) - t18;
t10 = -pkin(9) * t45 + (-pkin(3) - pkin(4)) * t35 + t70;
t31 = t35 * qJ(4);
t36 = t50 * t73 + t62 * t69;
t23 = -t37 * t58 + t59 * t41;
t68 = pkin(2) * t35 + t23;
t12 = t31 + (pkin(9) - qJ(3)) * t36 + t68;
t61 = sin(qJ(5));
t65 = cos(qJ(5));
t6 = t61 * t10 + t65 * t12;
t19 = -t62 * t34 + t37 * t71 + t41 * t72;
t17 = -t45 * pkin(2) + qJD(3) - t19;
t5 = t10 * t65 - t12 * t61;
t25 = t36 * t65 - t45 * t61;
t13 = t36 * pkin(3) - t45 * qJ(4) + t17;
t16 = -qJ(3) * t36 + t68;
t9 = -pkin(4) * t36 - t13;
t67 = V_base(3) ^ 2;
t64 = cos(qJ(6));
t60 = sin(qJ(6));
t33 = qJD(5) - t35;
t26 = t36 * t61 + t45 * t65;
t24 = qJD(6) - t25;
t22 = t26 * t64 + t33 * t60;
t21 = -t26 * t60 + t33 * t64;
t15 = -pkin(3) * t35 + t70;
t14 = t16 + t31;
t7 = -pkin(5) * t25 - pkin(10) * t26 + t9;
t4 = pkin(10) * t33 + t6;
t3 = -pkin(5) * t33 - t5;
t2 = t4 * t64 + t60 * t7;
t1 = -t4 * t60 + t64 * t7;
t8 = m(2) * (t46 ^ 2 + t47 ^ 2 + t67) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t67) / 0.2e1 + m(3) * (t19 ^ 2 + t20 ^ 2 + t23 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t46 * mrSges(2,1) - t47 * mrSges(2,2) + Ifges(2,3) * t57 / 0.2e1) * t57 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t33 / 0.2e1) * t33 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t46 * mrSges(2,3) + Ifges(2,5) * t57 + Ifges(2,1) * t50 / 0.2e1) * t50 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t33 + Ifges(6,1) * t26 / 0.2e1) * t26 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t22 / 0.2e1) * t22 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t47 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,6) * t57 + Ifges(2,2) * t49 / 0.2e1) * t49 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t33 + Ifges(6,2) * t25 / 0.2e1) * t25 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t22 + Ifges(7,6) * t24 + Ifges(7,2) * t21 / 0.2e1) * t21 + (t19 * mrSges(3,1) - t20 * mrSges(3,2) + t17 * mrSges(4,2) + t15 * mrSges(5,2) - t18 * mrSges(4,3) - t13 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t45) * t45 + (t17 * mrSges(4,1) + t13 * mrSges(5,1) + t23 * mrSges(3,2) - t14 * mrSges(5,2) - t19 * mrSges(3,3) - t16 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t36 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t45) * t36 + (t23 * mrSges(3,1) + t18 * mrSges(4,1) - t15 * mrSges(5,1) - t16 * mrSges(4,2) - t20 * mrSges(3,3) + t14 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t35 + (Ifges(5,4) + Ifges(4,5) - Ifges(3,6)) * t45 + (-Ifges(3,4) - Ifges(4,6) + Ifges(5,6)) * t36) * t35;
T  = t8;

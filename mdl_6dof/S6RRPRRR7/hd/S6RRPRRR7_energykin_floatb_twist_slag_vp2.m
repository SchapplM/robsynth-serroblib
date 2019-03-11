% Calculate kinetic energy for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:35
% EndTime: 2019-03-09 13:55:36
% DurationCPUTime: 1.10s
% Computational Cost: add. (2439->153), mult. (3057->213), div. (0->0), fcn. (2348->10), ass. (0->55)
t67 = sin(qJ(1));
t71 = cos(qJ(1));
t52 = t67 * V_base(5) + t71 * V_base(4);
t62 = V_base(6) + qJD(1);
t66 = sin(qJ(2));
t74 = cos(qJ(2));
t45 = t52 * t74 + t66 * t62;
t51 = -t67 * V_base(4) + t71 * V_base(5);
t50 = qJD(2) - t51;
t36 = -pkin(1) * t51 - pkin(7) * t52 + V_base(3);
t58 = V_base(5) * pkin(6) + V_base(1);
t59 = -V_base(4) * pkin(6) + V_base(2);
t47 = t71 * t58 + t67 * t59;
t41 = pkin(7) * t62 + t47;
t29 = t36 * t74 - t66 * t41;
t73 = qJD(3) - t29;
t19 = -t45 * pkin(8) + (-pkin(2) - pkin(3)) * t50 + t73;
t30 = t66 * t36 + t74 * t41;
t25 = t50 * qJ(3) + t30;
t44 = t52 * t66 - t62 * t74;
t22 = pkin(8) * t44 + t25;
t65 = sin(qJ(4));
t70 = cos(qJ(4));
t12 = t65 * t19 + t70 * t22;
t49 = qJD(4) - t50;
t10 = pkin(9) * t49 + t12;
t46 = -t67 * t58 + t71 * t59;
t40 = -t62 * pkin(1) - t46;
t26 = t44 * pkin(2) - t45 * qJ(3) + t40;
t23 = -pkin(3) * t44 - t26;
t33 = t44 * t70 - t65 * t45;
t34 = t44 * t65 + t45 * t70;
t15 = -pkin(4) * t33 - pkin(9) * t34 + t23;
t64 = sin(qJ(5));
t69 = cos(qJ(5));
t6 = t69 * t10 + t64 * t15;
t5 = -t10 * t64 + t69 * t15;
t11 = t19 * t70 - t65 * t22;
t32 = qJD(5) - t33;
t9 = -pkin(4) * t49 - t11;
t72 = V_base(3) ^ 2;
t68 = cos(qJ(6));
t63 = sin(qJ(6));
t31 = qJD(6) + t32;
t28 = t34 * t69 + t49 * t64;
t27 = -t34 * t64 + t49 * t69;
t24 = -t50 * pkin(2) + t73;
t17 = t27 * t63 + t28 * t68;
t16 = t27 * t68 - t28 * t63;
t7 = -pkin(5) * t27 + t9;
t4 = pkin(10) * t27 + t6;
t3 = pkin(5) * t32 - pkin(10) * t28 + t5;
t2 = t3 * t63 + t4 * t68;
t1 = t3 * t68 - t4 * t63;
t8 = (t46 * mrSges(2,1) - t47 * mrSges(2,2) + Ifges(2,3) * t62 / 0.2e1) * t62 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t32 + Ifges(6,2) * t27 / 0.2e1) * t27 + (t29 * mrSges(3,1) - t24 * mrSges(4,1) - t30 * mrSges(3,2) + t25 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t50) * t50 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t31 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t31 / 0.2e1) * t31 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t32 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t49 / 0.2e1) * t49 + (-V_base(3) * mrSges(2,1) + t47 * mrSges(2,3) + Ifges(2,4) * t52 + Ifges(2,6) * t62 + Ifges(2,2) * t51 / 0.2e1) * t51 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t46 * mrSges(2,3) + Ifges(2,5) * t62 + Ifges(2,1) * t52 / 0.2e1) * t52 + (t40 * mrSges(3,1) + t26 * mrSges(4,1) - t25 * mrSges(4,2) - t30 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t44 + (-Ifges(3,6) + Ifges(4,6)) * t50 + (-Ifges(3,4) + Ifges(4,5)) * t45) * t44 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t40 * mrSges(3,2) + t24 * mrSges(4,2) - t29 * mrSges(3,3) - t26 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t45 + (Ifges(4,4) + Ifges(3,5)) * t50) * t45 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t31 + Ifges(7,1) * t17 / 0.2e1) * t17 + m(2) * (t46 ^ 2 + t47 ^ 2 + t72) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t72) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t40 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t23 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (t23 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t49 + Ifges(5,1) * t34 / 0.2e1) * t34 + (-t23 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t34 + Ifges(5,6) * t49 + Ifges(5,2) * t33 / 0.2e1) * t33 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t32 / 0.2e1) * t32;
T  = t8;

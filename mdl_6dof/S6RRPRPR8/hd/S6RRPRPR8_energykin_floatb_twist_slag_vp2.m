% Calculate kinetic energy for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:26
% EndTime: 2019-03-09 10:49:27
% DurationCPUTime: 1.22s
% Computational Cost: add. (2829->153), mult. (3645->211), div. (0->0), fcn. (2876->10), ass. (0->56)
t73 = -pkin(4) - pkin(5);
t72 = cos(qJ(4));
t63 = sin(qJ(1));
t66 = cos(qJ(1));
t49 = -t63 * V_base(4) + t66 * V_base(5);
t50 = t63 * V_base(5) + t66 * V_base(4);
t36 = -pkin(1) * t49 - pkin(7) * t50 + V_base(3);
t55 = pkin(6) * V_base(5) + V_base(1);
t56 = -pkin(6) * V_base(4) + V_base(2);
t47 = t55 * t66 + t56 * t63;
t57 = V_base(6) + qJD(1);
t42 = pkin(7) * t57 + t47;
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t30 = t36 * t62 + t42 * t65;
t48 = qJD(2) - t49;
t27 = qJ(3) * t48 + t30;
t46 = -t55 * t63 + t56 * t66;
t41 = -pkin(1) * t57 - t46;
t44 = t50 * t62 - t57 * t65;
t45 = t50 * t65 + t57 * t62;
t28 = pkin(2) * t44 - qJ(3) * t45 + t41;
t58 = sin(pkin(10));
t59 = cos(pkin(10));
t18 = -t27 * t58 + t28 * t59;
t33 = t45 * t59 + t48 * t58;
t12 = pkin(3) * t44 - pkin(8) * t33 + t18;
t19 = t27 * t59 + t28 * t58;
t32 = -t45 * t58 + t48 * t59;
t15 = pkin(8) * t32 + t19;
t61 = sin(qJ(4));
t8 = t12 * t61 + t15 * t72;
t29 = t36 * t65 - t42 * t62;
t43 = qJD(4) + t44;
t6 = qJ(5) * t43 + t8;
t7 = t12 * t72 - t15 * t61;
t71 = pkin(2) * t48 - qJD(3) + t29;
t70 = qJD(5) - t7;
t69 = pkin(3) * t32 + t71;
t22 = t32 * t61 + t33 * t72;
t68 = qJ(5) * t22 + t69;
t67 = V_base(3) ^ 2;
t64 = cos(qJ(6));
t60 = sin(qJ(6));
t40 = qJD(6) - t43;
t21 = -t32 * t72 + t33 * t61;
t17 = t21 * t60 + t22 * t64;
t16 = t21 * t64 - t22 * t60;
t10 = pkin(4) * t21 - t68;
t9 = t21 * t73 + t68;
t5 = -pkin(4) * t43 + t70;
t4 = pkin(9) * t21 + t6;
t3 = -t22 * pkin(9) + t43 * t73 + t70;
t2 = t3 * t60 + t4 * t64;
t1 = t3 * t64 - t4 * t60;
t11 = m(2) * (t46 ^ 2 + t47 ^ 2 + t67) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t67) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t41 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t69 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t46 * mrSges(2,1) - t47 * mrSges(2,2) + Ifges(2,3) * t57 / 0.2e1) * t57 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t48 / 0.2e1) * t48 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t40 / 0.2e1) * t40 + (-t71 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,1) * t33 / 0.2e1) * t33 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t46 * mrSges(2,3) + Ifges(2,5) * t57 + Ifges(2,1) * t50 / 0.2e1) * t50 + (t41 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t48 + Ifges(3,1) * t45 / 0.2e1) * t45 + (t71 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t33 + Ifges(4,2) * t32 / 0.2e1) * t32 + (t9 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t40 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t47 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,6) * t57 + Ifges(2,2) * t49 / 0.2e1) * t49 + (-t9 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t40 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t7 * mrSges(5,1) - t5 * mrSges(6,1) - t8 * mrSges(5,2) + t6 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t43) * t43 + (-t69 * mrSges(5,2) + t5 * mrSges(6,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t22 + (Ifges(6,4) + Ifges(5,5)) * t43) * t22 + (t41 * mrSges(3,1) + t18 * mrSges(4,1) - t19 * mrSges(4,2) - t30 * mrSges(3,3) - Ifges(3,4) * t45 + Ifges(4,5) * t33 - Ifges(3,6) * t48 + Ifges(4,6) * t32 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t44) * t44 + (-t69 * mrSges(5,1) + t10 * mrSges(6,1) - t6 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t21 + (-Ifges(5,6) + Ifges(6,6)) * t43 + (-Ifges(5,4) + Ifges(6,5)) * t22) * t21;
T  = t11;

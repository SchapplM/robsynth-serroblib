% Calculate kinetic energy for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR14_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:27
% EndTime: 2019-03-09 11:32:28
% DurationCPUTime: 1.04s
% Computational Cost: add. (2709->157), mult. (4039->210), div. (0->0), fcn. (3276->10), ass. (0->60)
t74 = pkin(2) + pkin(9);
t73 = pkin(4) + pkin(10);
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t46 = t59 * V_base(5) + t61 * V_base(4);
t72 = pkin(8) * t46;
t71 = cos(qJ(2));
t70 = cos(qJ(4));
t58 = sin(qJ(2));
t45 = -t59 * V_base(4) + t61 * V_base(5);
t53 = V_base(6) + qJD(1);
t54 = sin(pkin(6));
t55 = cos(pkin(6));
t67 = t45 * t55 + t53 * t54;
t35 = t46 * t71 + t58 * t67;
t41 = -t45 * t54 + t53 * t55 + qJD(2);
t51 = V_base(5) * pkin(7) + V_base(1);
t52 = -V_base(4) * pkin(7) + V_base(2);
t43 = t61 * t51 + t59 * t52;
t33 = pkin(8) * t67 + t43;
t42 = -t51 * t59 + t61 * t52;
t36 = pkin(1) * t53 - t55 * t72 + t42;
t39 = -pkin(1) * t45 - t54 * t72 + V_base(3);
t68 = t55 * t71;
t69 = t54 * t71;
t20 = -t58 * t33 + t36 * t68 + t39 * t69;
t63 = qJD(3) - t20;
t12 = t35 * pkin(3) - t41 * t74 + t63;
t34 = -t45 * t68 + t46 * t58 - t53 * t69;
t24 = -t36 * t54 + t55 * t39;
t65 = -qJ(3) * t35 + t24;
t15 = t34 * t74 + t65;
t57 = sin(qJ(4));
t8 = t57 * t12 + t70 * t15;
t21 = t71 * t33 + (t36 * t55 + t39 * t54) * t58;
t19 = -t41 * qJ(3) - t21;
t32 = qJD(4) + t35;
t6 = -qJ(5) * t32 - t8;
t7 = t12 * t70 - t57 * t15;
t66 = qJD(5) - t7;
t16 = -pkin(3) * t34 - t19;
t27 = t57 * t34 + t41 * t70;
t64 = -qJ(5) * t27 + t16;
t62 = V_base(3) ^ 2;
t60 = cos(qJ(6));
t56 = sin(qJ(6));
t26 = -t34 * t70 + t41 * t57;
t25 = qJD(6) + t27;
t23 = t26 * t56 + t32 * t60;
t22 = t26 * t60 - t32 * t56;
t18 = -t41 * pkin(2) + t63;
t17 = pkin(2) * t34 + t65;
t10 = pkin(4) * t26 + t64;
t9 = t26 * t73 + t64;
t5 = -t32 * pkin(4) + t66;
t4 = -pkin(5) * t26 - t6;
t3 = t27 * pkin(5) - t32 * t73 + t66;
t2 = t3 * t56 + t60 * t9;
t1 = t3 * t60 - t56 * t9;
t11 = m(2) * (t42 ^ 2 + t43 ^ 2 + t62) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t62) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t24 ^ 2) / 0.2e1 + m(5) * (t16 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t42 * mrSges(2,1) - t43 * mrSges(2,2) + Ifges(2,3) * t53 / 0.2e1) * t53 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t42 * mrSges(2,3) + Ifges(2,5) * t53 + Ifges(2,1) * t46 / 0.2e1) * t46 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t43 * mrSges(2,3) + Ifges(2,4) * t46 + Ifges(2,6) * t53 + Ifges(2,2) * t45 / 0.2e1) * t45 + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t23 + Ifges(7,6) * t25 + Ifges(7,2) * t22 / 0.2e1) * t22 + (t20 * mrSges(3,1) - t21 * mrSges(3,2) + t18 * mrSges(4,2) - t19 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t41) * t41 + (t7 * mrSges(5,1) - t8 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t32) * t32 + (t18 * mrSges(4,1) + t24 * mrSges(3,2) - t20 * mrSges(3,3) - t17 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t35 + (-Ifges(4,4) + Ifges(3,5)) * t41) * t35 + (t5 * mrSges(6,1) + t16 * mrSges(5,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t27 + (-Ifges(6,4) + Ifges(5,5)) * t32) * t27 + (t24 * mrSges(3,1) + t19 * mrSges(4,1) - t17 * mrSges(4,2) - t21 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t34 + (Ifges(4,5) - Ifges(3,6)) * t41 + (-Ifges(3,4) - Ifges(4,6)) * t35) * t34 + (t16 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t26 + (Ifges(6,5) - Ifges(5,6)) * t32 + (-Ifges(5,4) - Ifges(6,6)) * t27) * t26;
T  = t11;

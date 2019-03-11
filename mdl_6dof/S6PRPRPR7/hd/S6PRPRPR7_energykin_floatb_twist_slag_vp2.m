% Calculate kinetic energy for
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:05
% EndTime: 2019-03-08 19:51:06
% DurationCPUTime: 0.98s
% Computational Cost: add. (2365->157), mult. (4039->209), div. (0->0), fcn. (3276->10), ass. (0->59)
t73 = pkin(2) + pkin(8);
t72 = pkin(4) + pkin(9);
t54 = sin(pkin(10));
t56 = cos(pkin(10));
t46 = t54 * V_base(5) + t56 * V_base(4);
t71 = pkin(7) * t46;
t70 = cos(qJ(2));
t69 = cos(qJ(4));
t60 = sin(qJ(2));
t45 = -t54 * V_base(4) + t56 * V_base(5);
t55 = sin(pkin(6));
t57 = cos(pkin(6));
t64 = t45 * t57 + t55 * V_base(6);
t35 = t46 * t70 + t60 * t64;
t41 = -t45 * t55 + t57 * V_base(6) + qJD(2);
t51 = V_base(5) * qJ(1) + V_base(1);
t52 = -V_base(4) * qJ(1) + V_base(2);
t43 = t56 * t51 + t54 * t52;
t32 = pkin(7) * t64 + t43;
t42 = -t51 * t54 + t56 * t52;
t36 = V_base(6) * pkin(1) - t57 * t71 + t42;
t53 = V_base(3) + qJD(1);
t39 = -pkin(1) * t45 - t55 * t71 + t53;
t67 = t57 * t70;
t68 = t55 * t70;
t20 = -t60 * t32 + t36 * t67 + t39 * t68;
t62 = qJD(3) - t20;
t12 = t35 * pkin(3) - t41 * t73 + t62;
t34 = -t45 * t67 + t46 * t60 - t68 * V_base(6);
t24 = -t36 * t55 + t57 * t39;
t65 = -qJ(3) * t35 + t24;
t15 = t34 * t73 + t65;
t59 = sin(qJ(4));
t8 = t59 * t12 + t69 * t15;
t21 = t70 * t32 + (t36 * t57 + t39 * t55) * t60;
t19 = -t41 * qJ(3) - t21;
t33 = qJD(4) + t35;
t6 = -qJ(5) * t33 - t8;
t7 = t12 * t69 - t59 * t15;
t66 = qJD(5) - t7;
t16 = -pkin(3) * t34 - t19;
t27 = t59 * t34 + t41 * t69;
t63 = -qJ(5) * t27 + t16;
t61 = cos(qJ(6));
t58 = sin(qJ(6));
t26 = -t34 * t69 + t41 * t59;
t25 = qJD(6) + t27;
t23 = t58 * t26 + t33 * t61;
t22 = t26 * t61 - t58 * t33;
t18 = -t41 * pkin(2) + t62;
t17 = pkin(2) * t34 + t65;
t10 = pkin(4) * t26 + t63;
t9 = t26 * t72 + t63;
t5 = -t33 * pkin(4) + t66;
t4 = -pkin(5) * t26 - t6;
t3 = t27 * pkin(5) - t33 * t72 + t66;
t2 = t58 * t3 + t61 * t9;
t1 = t3 * t61 - t58 * t9;
t11 = m(2) * (t42 ^ 2 + t43 ^ 2 + t53 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t24 ^ 2) / 0.2e1 + m(5) * (t16 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t53 * mrSges(2,2) - t42 * mrSges(2,3) + Ifges(2,1) * t46 / 0.2e1) * t46 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t53 * mrSges(2,1) + t43 * mrSges(2,3) + Ifges(2,4) * t46 + Ifges(2,2) * t45 / 0.2e1) * t45 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t23 / 0.2e1) * t23 + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t23 + Ifges(7,6) * t25 + Ifges(7,2) * t22 / 0.2e1) * t22 + (t20 * mrSges(3,1) - t21 * mrSges(3,2) + t18 * mrSges(4,2) - t19 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t41) * t41 + (t7 * mrSges(5,1) - t8 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t33) * t33 + (t18 * mrSges(4,1) + t24 * mrSges(3,2) - t20 * mrSges(3,3) - t17 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t35 + (-Ifges(4,4) + Ifges(3,5)) * t41) * t35 + (t5 * mrSges(6,1) + t16 * mrSges(5,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t27 + (-Ifges(6,4) + Ifges(5,5)) * t33) * t27 + (V_base(2) * mrSges(1,1) + t42 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t43 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t46 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t45 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t24 * mrSges(3,1) + t19 * mrSges(4,1) - t17 * mrSges(4,2) - t21 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t34 + (Ifges(4,5) - Ifges(3,6)) * t41 + (-Ifges(3,4) - Ifges(4,6)) * t35) * t34 + (t16 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t26 + (Ifges(6,5) - Ifges(5,6)) * t33 + (-Ifges(5,4) - Ifges(6,6)) * t27) * t26;
T  = t11;

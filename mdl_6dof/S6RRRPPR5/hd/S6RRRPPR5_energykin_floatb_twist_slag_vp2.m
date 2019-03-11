% Calculate kinetic energy for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:09
% EndTime: 2019-03-09 15:38:11
% DurationCPUTime: 1.44s
% Computational Cost: add. (6471->161), mult. (9743->236), div. (0->0), fcn. (8340->14), ass. (0->65)
t62 = pkin(7) * V_base(5) + V_base(1);
t63 = -pkin(7) * V_base(4) + V_base(2);
t73 = sin(qJ(1));
t77 = cos(qJ(1));
t55 = -t62 * t73 + t63 * t77;
t64 = V_base(6) + qJD(1);
t69 = cos(pkin(6));
t58 = t73 * V_base(5) + t77 * V_base(4);
t85 = pkin(8) * t58;
t50 = pkin(1) * t64 - t69 * t85 + t55;
t57 = -t73 * V_base(4) + t77 * V_base(5);
t67 = sin(pkin(6));
t53 = -pkin(1) * t57 - t67 * t85 + V_base(3);
t86 = t50 * t69 + t53 * t67;
t56 = t62 * t77 + t63 * t73;
t79 = t57 * t69 + t64 * t67;
t47 = pkin(8) * t79 + t56;
t72 = sin(qJ(2));
t76 = cos(qJ(2));
t37 = -t47 * t72 + t76 * t86;
t48 = -t58 * t72 + t76 * t79;
t39 = -t50 * t67 + t53 * t69;
t49 = t58 * t76 + t72 * t79;
t30 = -pkin(2) * t48 - pkin(9) * t49 + t39;
t38 = t76 * t47 + t72 * t86;
t54 = -t57 * t67 + t64 * t69 + qJD(2);
t33 = pkin(9) * t54 + t38;
t71 = sin(qJ(3));
t75 = cos(qJ(3));
t23 = t30 * t75 - t33 * t71;
t42 = t49 * t75 + t54 * t71;
t46 = qJD(3) - t48;
t17 = pkin(3) * t46 - qJ(4) * t42 + t23;
t24 = t30 * t71 + t33 * t75;
t41 = -t49 * t71 + t54 * t75;
t22 = qJ(4) * t41 + t24;
t66 = sin(pkin(11));
t81 = cos(pkin(11));
t12 = t17 * t66 + t22 * t81;
t10 = qJ(5) * t46 + t12;
t32 = -pkin(2) * t54 - t37;
t25 = -pkin(3) * t41 + qJD(4) + t32;
t35 = -t41 * t81 + t42 * t66;
t36 = t41 * t66 + t42 * t81;
t15 = pkin(4) * t35 - qJ(5) * t36 + t25;
t65 = sin(pkin(12));
t68 = cos(pkin(12));
t6 = t10 * t68 + t15 * t65;
t5 = -t10 * t65 + t15 * t68;
t11 = t17 * t81 - t22 * t66;
t9 = -pkin(4) * t46 + qJD(5) - t11;
t78 = V_base(3) ^ 2;
t74 = cos(qJ(6));
t70 = sin(qJ(6));
t34 = qJD(6) + t35;
t27 = t36 * t68 + t46 * t65;
t26 = -t36 * t65 + t46 * t68;
t19 = t26 * t70 + t27 * t74;
t18 = t26 * t74 - t27 * t70;
t7 = -pkin(5) * t26 + t9;
t4 = pkin(10) * t26 + t6;
t3 = pkin(5) * t35 - pkin(10) * t27 + t5;
t2 = t3 * t70 + t4 * t74;
t1 = t3 * t74 - t4 * t70;
t8 = (t39 * mrSges(3,2) - t37 * mrSges(3,3) + Ifges(3,5) * t54 + Ifges(3,1) * t49 / 0.2e1) * t49 + (-V_base(3) * mrSges(2,1) + t56 * mrSges(2,3) + Ifges(2,4) * t58 + Ifges(2,6) * t64 + Ifges(2,2) * t57 / 0.2e1) * t57 + (t25 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t36 + Ifges(6,5) * t27 - Ifges(5,6) * t46 + Ifges(6,6) * t26 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t35) * t35 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,2) * t26 / 0.2e1) * t26 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t39 * mrSges(3,1) + t38 * mrSges(3,3) + Ifges(3,4) * t49 + Ifges(3,6) * t54 + Ifges(3,2) * t48 / 0.2e1) * t48 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(2,2) - t55 * mrSges(2,3) + Ifges(2,5) * t64 + Ifges(2,1) * t58 / 0.2e1) * t58 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t27 / 0.2e1) * t27 + (-t32 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t42 + Ifges(4,2) * t41 / 0.2e1) * t41 + (t55 * mrSges(2,1) - t56 * mrSges(2,2) + Ifges(2,3) * t64 / 0.2e1) * t64 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t34 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t37 * mrSges(3,1) - t38 * mrSges(3,2) + Ifges(3,3) * t54 / 0.2e1) * t54 + (t32 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,1) * t42 / 0.2e1) * t42 + (t23 * mrSges(4,1) + t11 * mrSges(5,1) - t24 * mrSges(4,2) - t12 * mrSges(5,2) + Ifges(4,5) * t42 + Ifges(5,5) * t36 + Ifges(4,6) * t41 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t46) * t46 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (t25 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,1) * t36 / 0.2e1) * t36 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t37 ^ 2 + t38 ^ 2 + t39 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t78) / 0.2e1 + m(2) * (t55 ^ 2 + t56 ^ 2 + t78) / 0.2e1;
T  = t8;

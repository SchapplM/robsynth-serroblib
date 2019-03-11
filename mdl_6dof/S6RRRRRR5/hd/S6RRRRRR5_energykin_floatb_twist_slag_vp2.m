% Calculate kinetic energy for
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:22
% EndTime: 2019-03-10 03:56:23
% DurationCPUTime: 1.69s
% Computational Cost: add. (6975->161), mult. (10451->240), div. (0->0), fcn. (8960->14), ass. (0->67)
t63 = V_base(5) * pkin(7) + V_base(1);
t64 = -V_base(4) * pkin(7) + V_base(2);
t73 = sin(qJ(1));
t79 = cos(qJ(1));
t55 = -t63 * t73 + t79 * t64;
t65 = V_base(6) + qJD(1);
t67 = cos(pkin(6));
t59 = t73 * V_base(5) + t79 * V_base(4);
t85 = pkin(8) * t59;
t50 = pkin(1) * t65 - t67 * t85 + t55;
t58 = -t73 * V_base(4) + t79 * V_base(5);
t66 = sin(pkin(6));
t53 = -pkin(1) * t58 - t66 * t85 + V_base(3);
t86 = t50 * t67 + t53 * t66;
t56 = t79 * t63 + t73 * t64;
t81 = t58 * t67 + t65 * t66;
t47 = pkin(8) * t81 + t56;
t72 = sin(qJ(2));
t78 = cos(qJ(2));
t36 = -t72 * t47 + t78 * t86;
t48 = -t72 * t59 + t78 * t81;
t38 = -t50 * t66 + t67 * t53;
t49 = t59 * t78 + t72 * t81;
t30 = -pkin(2) * t48 - pkin(9) * t49 + t38;
t37 = t78 * t47 + t86 * t72;
t54 = -t58 * t66 + t65 * t67 + qJD(2);
t33 = pkin(9) * t54 + t37;
t71 = sin(qJ(3));
t77 = cos(qJ(3));
t22 = t77 * t30 - t33 * t71;
t40 = t49 * t77 + t54 * t71;
t46 = qJD(3) - t48;
t16 = pkin(3) * t46 - pkin(10) * t40 + t22;
t23 = t71 * t30 + t77 * t33;
t39 = -t49 * t71 + t54 * t77;
t18 = pkin(10) * t39 + t23;
t70 = sin(qJ(4));
t76 = cos(qJ(4));
t13 = t70 * t16 + t76 * t18;
t34 = t39 * t76 - t40 * t70;
t11 = pkin(11) * t34 + t13;
t69 = sin(qJ(5));
t75 = cos(qJ(5));
t12 = t76 * t16 - t18 * t70;
t35 = t39 * t70 + t40 * t76;
t45 = qJD(4) + t46;
t8 = pkin(4) * t45 - pkin(11) * t35 + t12;
t6 = t75 * t11 + t69 * t8;
t5 = -t11 * t69 + t75 * t8;
t25 = t34 * t75 - t35 * t69;
t32 = -pkin(2) * t54 - t36;
t27 = -pkin(3) * t39 + t32;
t19 = -pkin(4) * t34 + t27;
t80 = V_base(3) ^ 2;
t74 = cos(qJ(6));
t68 = sin(qJ(6));
t41 = qJD(5) + t45;
t26 = t34 * t69 + t35 * t75;
t24 = qJD(6) - t25;
t21 = t26 * t74 + t41 * t68;
t20 = -t26 * t68 + t41 * t74;
t10 = -pkin(5) * t25 - pkin(12) * t26 + t19;
t4 = pkin(12) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t10 * t68 + t4 * t74;
t1 = t10 * t74 - t4 * t68;
t7 = (-V_base(3) * mrSges(2,1) + t56 * mrSges(2,3) + Ifges(2,4) * t59 + Ifges(2,6) * t65 + Ifges(2,2) * t58 / 0.2e1) * t58 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t41 / 0.2e1) * t41 + (t36 * mrSges(3,1) - t37 * mrSges(3,2) + Ifges(3,3) * t54 / 0.2e1) * t54 + (-t27 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t35 + Ifges(5,6) * t45 + Ifges(5,2) * t34 / 0.2e1) * t34 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t21 / 0.2e1) * t21 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t24 + Ifges(7,2) * t20 / 0.2e1) * t20 + (-t19 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t41 + Ifges(6,2) * t25 / 0.2e1) * t25 + (t19 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t41 + Ifges(6,1) * t26 / 0.2e1) * t26 + (-t38 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t49 + Ifges(3,6) * t54 + Ifges(3,2) * t48 / 0.2e1) * t48 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24 + (V_base(3) * mrSges(2,2) - t55 * mrSges(2,3) + Ifges(2,5) * t65 + Ifges(2,1) * t59 / 0.2e1) * t59 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t27 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t45 + Ifges(5,1) * t35 / 0.2e1) * t35 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t32 * mrSges(4,2) - t22 * mrSges(4,3) + Ifges(4,5) * t46 + Ifges(4,1) * t40 / 0.2e1) * t40 + (t22 * mrSges(4,1) - t23 * mrSges(4,2) + Ifges(4,3) * t46 / 0.2e1) * t46 + (t38 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,5) * t54 + Ifges(3,1) * t49 / 0.2e1) * t49 + (t55 * mrSges(2,1) - t56 * mrSges(2,2) + Ifges(2,3) * t65 / 0.2e1) * t65 + m(2) * (t55 ^ 2 + t56 ^ 2 + t80) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t80) / 0.2e1 + (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t45 / 0.2e1) * t45 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t22 ^ 2 + t23 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + (-t32 * mrSges(4,1) + t23 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,6) * t46 + Ifges(4,2) * t39 / 0.2e1) * t39;
T  = t7;

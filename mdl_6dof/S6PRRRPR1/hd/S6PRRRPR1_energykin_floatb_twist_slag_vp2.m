% Calculate kinetic energy for
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:04
% EndTime: 2019-03-08 23:00:05
% DurationCPUTime: 1.51s
% Computational Cost: add. (6041->161), mult. (10451->237), div. (0->0), fcn. (8960->14), ass. (0->65)
t62 = V_base(5) * qJ(1) + V_base(1);
t63 = -V_base(4) * qJ(1) + V_base(2);
t66 = sin(pkin(11));
t69 = cos(pkin(11));
t54 = -t62 * t66 + t69 * t63;
t70 = cos(pkin(6));
t58 = t66 * V_base(5) + t69 * V_base(4);
t83 = pkin(7) * t58;
t49 = V_base(6) * pkin(1) - t70 * t83 + t54;
t57 = -t66 * V_base(4) + t69 * V_base(5);
t64 = V_base(3) + qJD(1);
t67 = sin(pkin(6));
t52 = -pkin(1) * t57 - t67 * t83 + t64;
t84 = t49 * t70 + t52 * t67;
t55 = t69 * t62 + t66 * t63;
t79 = t57 * t70 + t67 * V_base(6);
t45 = pkin(7) * t79 + t55;
t74 = sin(qJ(2));
t78 = cos(qJ(2));
t36 = -t74 * t45 + t78 * t84;
t47 = -t74 * t58 + t78 * t79;
t38 = -t49 * t67 + t70 * t52;
t48 = t58 * t78 + t74 * t79;
t30 = -pkin(2) * t47 - pkin(8) * t48 + t38;
t37 = t78 * t45 + t84 * t74;
t53 = -t57 * t67 + t70 * V_base(6) + qJD(2);
t33 = pkin(8) * t53 + t37;
t73 = sin(qJ(3));
t77 = cos(qJ(3));
t22 = t77 * t30 - t33 * t73;
t40 = t48 * t77 + t53 * t73;
t46 = qJD(3) - t47;
t16 = pkin(3) * t46 - pkin(9) * t40 + t22;
t23 = t73 * t30 + t77 * t33;
t39 = -t48 * t73 + t53 * t77;
t18 = pkin(9) * t39 + t23;
t72 = sin(qJ(4));
t76 = cos(qJ(4));
t13 = t72 * t16 + t76 * t18;
t34 = t39 * t76 - t40 * t72;
t10 = qJ(5) * t34 + t13;
t65 = sin(pkin(12));
t68 = cos(pkin(12));
t12 = t76 * t16 - t18 * t72;
t35 = t39 * t72 + t40 * t76;
t44 = qJD(4) + t46;
t8 = pkin(4) * t44 - qJ(5) * t35 + t12;
t6 = t68 * t10 + t65 * t8;
t5 = -t10 * t65 + t68 * t8;
t25 = t34 * t68 - t35 * t65;
t32 = -t53 * pkin(2) - t36;
t27 = -t39 * pkin(3) + t32;
t19 = -t34 * pkin(4) + qJD(5) + t27;
t75 = cos(qJ(6));
t71 = sin(qJ(6));
t26 = t34 * t65 + t35 * t68;
t24 = qJD(6) - t25;
t21 = t26 * t75 + t44 * t71;
t20 = -t26 * t71 + t44 * t75;
t11 = -t25 * pkin(5) - t26 * pkin(10) + t19;
t4 = pkin(10) * t44 + t6;
t3 = -pkin(5) * t44 - t5;
t2 = t11 * t71 + t4 * t75;
t1 = t11 * t75 - t4 * t71;
t7 = (t19 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t26 / 0.2e1) * t26 + (-t38 * mrSges(3,1) + t37 * mrSges(3,3) + Ifges(3,4) * t48 + Ifges(3,6) * t53 + Ifges(3,2) * t47 / 0.2e1) * t47 + (t36 * mrSges(3,1) - t37 * mrSges(3,2) + Ifges(3,3) * t53 / 0.2e1) * t53 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t13 * mrSges(5,2) - t6 * mrSges(6,2) + Ifges(5,5) * t35 + Ifges(6,5) * t26 + Ifges(5,6) * t34 + Ifges(6,6) * t25 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t44) * t44 + (-t32 * mrSges(4,1) + t23 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,6) * t46 + Ifges(4,2) * t39 / 0.2e1) * t39 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t24 + Ifges(7,2) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t32 * mrSges(4,2) - t22 * mrSges(4,3) + Ifges(4,5) * t46 + Ifges(4,1) * t40 / 0.2e1) * t40 + (-t27 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t35 + Ifges(5,2) * t34 / 0.2e1) * t34 + (t22 * mrSges(4,1) - t23 * mrSges(4,2) + Ifges(4,3) * t46 / 0.2e1) * t46 + (-t19 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,2) * t25 / 0.2e1) * t25 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t21 / 0.2e1) * t21 + (t27 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t35 / 0.2e1) * t35 + (V_base(2) * mrSges(1,1) + t54 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t55 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t58 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t57 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-t64 * mrSges(2,1) + t55 * mrSges(2,3) + Ifges(2,4) * t58 + Ifges(2,2) * t57 / 0.2e1) * t57 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t54 ^ 2 + t55 ^ 2 + t64 ^ 2) / 0.2e1 + (t38 * mrSges(3,2) - t36 * mrSges(3,3) + t53 * Ifges(3,5) + Ifges(3,1) * t48 / 0.2e1) * t48 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (t64 * mrSges(2,2) - t54 * mrSges(2,3) + Ifges(2,1) * t58 / 0.2e1) * t58 + m(6) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t22 ^ 2 + t23 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24;
T  = t7;

% Calculate kinetic energy for
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:32
% EndTime: 2019-03-08 19:24:34
% DurationCPUTime: 1.38s
% Computational Cost: add. (5459->161), mult. (9889->235), div. (0->0), fcn. (8486->14), ass. (0->64)
t61 = V_base(5) * qJ(1) + V_base(1);
t62 = -V_base(4) * qJ(1) + V_base(2);
t66 = sin(pkin(10));
t70 = cos(pkin(10));
t54 = -t61 * t66 + t70 * t62;
t71 = cos(pkin(6));
t57 = t66 * V_base(5) + t70 * V_base(4);
t81 = pkin(7) * t57;
t48 = V_base(6) * pkin(1) - t71 * t81 + t54;
t56 = -t66 * V_base(4) + t70 * V_base(5);
t63 = V_base(3) + qJD(1);
t67 = sin(pkin(6));
t52 = -pkin(1) * t56 - t67 * t81 + t63;
t82 = t48 * t71 + t52 * t67;
t55 = t70 * t61 + t66 * t62;
t78 = t56 * t71 + t67 * V_base(6);
t45 = t78 * pkin(7) + t55;
t74 = sin(qJ(2));
t77 = cos(qJ(2));
t33 = -t45 * t74 + t82 * t77;
t47 = t57 * t77 + t78 * t74;
t53 = -t56 * t67 + t71 * V_base(6) + qJD(2);
t29 = pkin(2) * t53 - qJ(3) * t47 + t33;
t34 = t77 * t45 + t82 * t74;
t46 = -t74 * t57 + t78 * t77;
t32 = qJ(3) * t46 + t34;
t65 = sin(pkin(11));
t69 = cos(pkin(11));
t19 = t65 * t29 + t69 * t32;
t17 = pkin(8) * t53 + t19;
t41 = -t48 * t67 + t71 * t52;
t35 = -pkin(2) * t46 + qJD(3) + t41;
t39 = t46 * t69 - t47 * t65;
t40 = t46 * t65 + t47 * t69;
t24 = -pkin(3) * t39 - pkin(8) * t40 + t35;
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t13 = t76 * t17 + t73 * t24;
t36 = -t40 * t73 + t53 * t76;
t11 = qJ(5) * t36 + t13;
t64 = sin(pkin(12));
t68 = cos(pkin(12));
t12 = -t17 * t73 + t76 * t24;
t37 = t40 * t76 + t53 * t73;
t38 = qJD(4) - t39;
t9 = pkin(4) * t38 - qJ(5) * t37 + t12;
t6 = t68 * t11 + t64 * t9;
t18 = t29 * t69 - t65 * t32;
t5 = -t11 * t64 + t68 * t9;
t26 = t36 * t68 - t37 * t64;
t16 = -pkin(3) * t53 - t18;
t14 = -pkin(4) * t36 + qJD(5) + t16;
t75 = cos(qJ(6));
t72 = sin(qJ(6));
t27 = t36 * t64 + t37 * t68;
t25 = qJD(6) - t26;
t21 = t27 * t75 + t38 * t72;
t20 = -t27 * t72 + t38 * t75;
t7 = -pkin(5) * t26 - pkin(9) * t27 + t14;
t4 = pkin(9) * t38 + t6;
t3 = -pkin(5) * t38 - t5;
t2 = t4 * t75 + t7 * t72;
t1 = -t4 * t72 + t7 * t75;
t8 = (t33 * mrSges(3,1) + t18 * mrSges(4,1) - t34 * mrSges(3,2) - t19 * mrSges(4,2) + Ifges(3,5) * t47 + Ifges(4,5) * t40 + Ifges(3,6) * t46 + Ifges(4,6) * t39 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t53) * t53 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (-t35 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,2) * t39 / 0.2e1) * t39 + (V_base(2) * mrSges(1,1) + t54 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t55 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t57 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t56 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-t41 * mrSges(3,1) + t34 * mrSges(3,3) + Ifges(3,4) * t47 + Ifges(3,2) * t46 / 0.2e1) * t46 + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t13 * mrSges(5,2) - t6 * mrSges(6,2) + Ifges(5,5) * t37 + Ifges(6,5) * t27 + Ifges(5,6) * t36 + Ifges(6,6) * t26 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t38) * t38 + (t16 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t37 / 0.2e1) * t37 + (-t16 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,2) * t36 / 0.2e1) * t36 + (-t63 * mrSges(2,1) + t55 * mrSges(2,3) + Ifges(2,4) * t57 + Ifges(2,2) * t56 / 0.2e1) * t56 + (t35 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,1) * t40 / 0.2e1) * t40 + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t63 * mrSges(2,2) - t54 * mrSges(2,3) + Ifges(2,1) * t57 / 0.2e1) * t57 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t21 / 0.2e1) * t21 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t25 + Ifges(7,2) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t41 ^ 2) / 0.2e1 + m(2) * (t54 ^ 2 + t55 ^ 2 + t63 ^ 2) / 0.2e1 + (t41 * mrSges(3,2) - t33 * mrSges(3,3) + Ifges(3,1) * t47 / 0.2e1) * t47;
T  = t8;

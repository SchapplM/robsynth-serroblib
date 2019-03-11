% Calculate kinetic energy for
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:13:49
% EndTime: 2019-03-08 19:13:50
% DurationCPUTime: 1.36s
% Computational Cost: add. (5445->161), mult. (9889->235), div. (0->0), fcn. (8486->14), ass. (0->64)
t62 = V_base(5) * qJ(1) + V_base(1);
t63 = -V_base(4) * qJ(1) + V_base(2);
t67 = sin(pkin(10));
t71 = cos(pkin(10));
t55 = -t62 * t67 + t71 * t63;
t72 = cos(pkin(6));
t58 = t67 * V_base(5) + t71 * V_base(4);
t82 = pkin(7) * t58;
t49 = V_base(6) * pkin(1) - t72 * t82 + t55;
t57 = -t67 * V_base(4) + t71 * V_base(5);
t64 = V_base(3) + qJD(1);
t68 = sin(pkin(6));
t53 = -pkin(1) * t57 - t68 * t82 + t64;
t83 = t49 * t72 + t53 * t68;
t56 = t71 * t62 + t67 * t63;
t79 = t57 * t72 + t68 * V_base(6);
t46 = t79 * pkin(7) + t56;
t75 = sin(qJ(2));
t78 = cos(qJ(2));
t33 = -t46 * t75 + t83 * t78;
t48 = t58 * t78 + t79 * t75;
t54 = -t57 * t68 + t72 * V_base(6) + qJD(2);
t29 = pkin(2) * t54 - qJ(3) * t48 + t33;
t34 = t78 * t46 + t83 * t75;
t47 = -t75 * t58 + t79 * t78;
t32 = qJ(3) * t47 + t34;
t66 = sin(pkin(11));
t70 = cos(pkin(11));
t19 = t66 * t29 + t70 * t32;
t17 = qJ(4) * t54 + t19;
t41 = -t49 * t68 + t72 * t53;
t35 = -pkin(2) * t47 + qJD(3) + t41;
t39 = -t70 * t47 + t48 * t66;
t40 = t47 * t66 + t48 * t70;
t24 = pkin(3) * t39 - qJ(4) * t40 + t35;
t65 = sin(pkin(12));
t69 = cos(pkin(12));
t13 = t69 * t17 + t65 * t24;
t36 = -t40 * t65 + t54 * t69;
t11 = pkin(8) * t36 + t13;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t12 = -t17 * t65 + t69 * t24;
t37 = t40 * t69 + t54 * t65;
t9 = pkin(4) * t39 - pkin(8) * t37 + t12;
t6 = t77 * t11 + t74 * t9;
t18 = t29 * t70 - t66 * t32;
t5 = -t11 * t74 + t77 * t9;
t26 = t36 * t77 - t37 * t74;
t16 = -pkin(3) * t54 + qJD(4) - t18;
t14 = -pkin(4) * t36 + t16;
t76 = cos(qJ(6));
t73 = sin(qJ(6));
t38 = qJD(5) + t39;
t27 = t36 * t74 + t37 * t77;
t25 = qJD(6) - t26;
t21 = t27 * t76 + t38 * t73;
t20 = -t27 * t73 + t38 * t76;
t7 = -pkin(5) * t26 - pkin(9) * t27 + t14;
t4 = pkin(9) * t38 + t6;
t3 = -pkin(5) * t38 - t5;
t2 = t4 * t76 + t7 * t73;
t1 = -t4 * t73 + t7 * t76;
t8 = (V_base(2) * mrSges(1,1) + t55 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t56 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t58 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t57 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (t35 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,1) * t40 / 0.2e1) * t40 + (t33 * mrSges(3,1) + t18 * mrSges(4,1) - t34 * mrSges(3,2) - t19 * mrSges(4,2) + Ifges(3,5) * t48 + Ifges(4,5) * t40 + Ifges(3,6) * t47 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t54) * t54 + (t35 * mrSges(4,1) + t12 * mrSges(5,1) - t13 * mrSges(5,2) - t19 * mrSges(4,3) - Ifges(4,4) * t40 + Ifges(5,5) * t37 - Ifges(4,6) * t54 + Ifges(5,6) * t36 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t39) * t39 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t41 * mrSges(3,1) + t34 * mrSges(3,3) + Ifges(3,4) * t48 + Ifges(3,2) * t47 / 0.2e1) * t47 + (-t16 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,2) * t36 / 0.2e1) * t36 + (t64 * mrSges(2,2) - t55 * mrSges(2,3) + Ifges(2,1) * t58 / 0.2e1) * t58 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t25 + Ifges(7,1) * t21 / 0.2e1) * t21 + (-t14 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,6) * t38 + Ifges(6,2) * t26 / 0.2e1) * t26 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t25 + Ifges(7,2) * t20 / 0.2e1) * t20 + (t41 * mrSges(3,2) - t33 * mrSges(3,3) + Ifges(3,1) * t48 / 0.2e1) * t48 + (-t64 * mrSges(2,1) + t56 * mrSges(2,3) + Ifges(2,4) * t58 + Ifges(2,2) * t57 / 0.2e1) * t57 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t38 / 0.2e1) * t38 + (t16 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t37 / 0.2e1) * t37 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t16 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t41 ^ 2) / 0.2e1 + (t14 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t38 + Ifges(6,1) * t27 / 0.2e1) * t27 + m(2) * (t55 ^ 2 + t56 ^ 2 + t64 ^ 2) / 0.2e1;
T  = t8;

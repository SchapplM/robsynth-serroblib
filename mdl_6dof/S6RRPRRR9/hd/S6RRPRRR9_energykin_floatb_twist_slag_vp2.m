% Calculate kinetic energy for
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:13
% EndTime: 2019-03-09 14:08:14
% DurationCPUTime: 1.58s
% Computational Cost: add. (6915->161), mult. (10451->238), div. (0->0), fcn. (8960->14), ass. (0->66)
t74 = sin(qJ(1));
t79 = cos(qJ(1));
t57 = -t74 * V_base(4) + t79 * V_base(5);
t65 = V_base(6) + qJD(1);
t67 = sin(pkin(6));
t69 = cos(pkin(6));
t81 = t57 * t69 + t65 * t67;
t63 = V_base(5) * pkin(7) + V_base(1);
t64 = -V_base(4) * pkin(7) + V_base(2);
t54 = -t63 * t74 + t79 * t64;
t58 = t74 * V_base(5) + t79 * V_base(4);
t87 = pkin(8) * t58;
t49 = pkin(1) * t65 - t69 * t87 + t54;
t52 = -pkin(1) * t57 - t67 * t87 + V_base(3);
t88 = t49 * t69 + t52 * t67;
t55 = t79 * t63 + t74 * t64;
t46 = t81 * pkin(8) + t55;
t73 = sin(qJ(2));
t78 = cos(qJ(2));
t36 = -t73 * t46 + t88 * t78;
t38 = -t49 * t67 + t69 * t52;
t47 = t58 * t73 - t81 * t78;
t48 = t58 * t78 + t81 * t73;
t30 = pkin(2) * t47 - qJ(3) * t48 + t38;
t37 = t78 * t46 + t88 * t73;
t53 = -t57 * t67 + t65 * t69 + qJD(2);
t33 = qJ(3) * t53 + t37;
t66 = sin(pkin(12));
t68 = cos(pkin(12));
t22 = t68 * t30 - t33 * t66;
t40 = t48 * t68 + t53 * t66;
t16 = pkin(3) * t47 - pkin(9) * t40 + t22;
t23 = t66 * t30 + t68 * t33;
t39 = -t48 * t66 + t53 * t68;
t18 = pkin(9) * t39 + t23;
t72 = sin(qJ(4));
t77 = cos(qJ(4));
t13 = t72 * t16 + t77 * t18;
t34 = t39 * t77 - t40 * t72;
t11 = pkin(10) * t34 + t13;
t71 = sin(qJ(5));
t76 = cos(qJ(5));
t12 = t77 * t16 - t18 * t72;
t35 = t39 * t72 + t40 * t77;
t45 = qJD(4) + t47;
t8 = pkin(4) * t45 - pkin(10) * t35 + t12;
t6 = t76 * t11 + t71 * t8;
t5 = -t11 * t71 + t76 * t8;
t25 = t34 * t76 - t35 * t71;
t32 = -pkin(2) * t53 + qJD(3) - t36;
t27 = -pkin(3) * t39 + t32;
t19 = -pkin(4) * t34 + t27;
t80 = V_base(3) ^ 2;
t75 = cos(qJ(6));
t70 = sin(qJ(6));
t44 = qJD(5) + t45;
t26 = t34 * t71 + t35 * t76;
t24 = qJD(6) - t25;
t21 = t26 * t75 + t44 * t70;
t20 = -t26 * t70 + t44 * t75;
t9 = -pkin(5) * t25 - pkin(11) * t26 + t19;
t4 = pkin(11) * t44 + t6;
t3 = -pkin(5) * t44 - t5;
t2 = t4 * t75 + t70 * t9;
t1 = -t4 * t70 + t75 * t9;
t7 = (t12 * mrSges(5,1) - t13 * mrSges(5,2) + Ifges(5,3) * t45 / 0.2e1) * t45 + (-t32 * mrSges(4,1) + t23 * mrSges(4,3) + Ifges(4,4) * t40 + Ifges(4,2) * t39 / 0.2e1) * t39 + (-t27 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t35 + Ifges(5,6) * t45 + Ifges(5,2) * t34 / 0.2e1) * t34 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t24 + Ifges(7,2) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(2,1) + t55 * mrSges(2,3) + Ifges(2,4) * t58 + Ifges(2,6) * t65 + Ifges(2,2) * t57 / 0.2e1) * t57 + (-t19 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t44 + Ifges(6,2) * t25 / 0.2e1) * t25 + (t19 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t44 + Ifges(6,1) * t26 / 0.2e1) * t26 + (t38 * mrSges(3,1) + t22 * mrSges(4,1) - t23 * mrSges(4,2) - t37 * mrSges(3,3) - Ifges(3,4) * t48 + Ifges(4,5) * t40 - Ifges(3,6) * t53 + Ifges(4,6) * t39 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t47) * t47 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t24 + Ifges(7,1) * t21 / 0.2e1) * t21 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t44 / 0.2e1) * t44 + (V_base(3) * mrSges(2,2) - t54 * mrSges(2,3) + Ifges(2,5) * t65 + Ifges(2,1) * t58 / 0.2e1) * t58 + (t27 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,5) * t45 + Ifges(5,1) * t35 / 0.2e1) * t35 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t54 * mrSges(2,1) - t55 * mrSges(2,2) + Ifges(2,3) * t65 / 0.2e1) * t65 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t24 / 0.2e1) * t24 + (t38 * mrSges(3,2) - t36 * mrSges(3,3) + Ifges(3,5) * t53 + Ifges(3,1) * t48 / 0.2e1) * t48 + (t36 * mrSges(3,1) - t37 * mrSges(3,2) + Ifges(3,3) * t53 / 0.2e1) * t53 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (t32 * mrSges(4,2) - t22 * mrSges(4,3) + Ifges(4,1) * t40 / 0.2e1) * t40 + m(5) * (t12 ^ 2 + t13 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t22 ^ 2 + t23 ^ 2 + t32 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t80) / 0.2e1 + m(2) * (t54 ^ 2 + t55 ^ 2 + t80) / 0.2e1;
T  = t7;

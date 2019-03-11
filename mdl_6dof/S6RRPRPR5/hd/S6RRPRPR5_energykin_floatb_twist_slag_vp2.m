% Calculate kinetic energy for
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:28:41
% EndTime: 2019-03-09 10:28:42
% DurationCPUTime: 1.44s
% Computational Cost: add. (6313->161), mult. (9957->236), div. (0->0), fcn. (8538->14), ass. (0->65)
t63 = V_base(5) * pkin(7) + V_base(1);
t64 = -V_base(4) * pkin(7) + V_base(2);
t75 = sin(qJ(1));
t78 = cos(qJ(1));
t56 = -t63 * t75 + t78 * t64;
t65 = V_base(6) + qJD(1);
t71 = cos(pkin(6));
t59 = t75 * V_base(5) + t78 * V_base(4);
t84 = pkin(8) * t59;
t49 = pkin(1) * t65 - t71 * t84 + t56;
t58 = -t75 * V_base(4) + t78 * V_base(5);
t68 = sin(pkin(6));
t53 = -pkin(1) * t58 - t68 * t84 + V_base(3);
t85 = t49 * t71 + t53 * t68;
t57 = t78 * t63 + t75 * t64;
t80 = t58 * t71 + t65 * t68;
t46 = pkin(8) * t80 + t57;
t74 = sin(qJ(2));
t77 = cos(qJ(2));
t33 = -t46 * t74 + t85 * t77;
t48 = t59 * t77 + t74 * t80;
t55 = -t58 * t68 + t65 * t71 + qJD(2);
t29 = pkin(2) * t55 - qJ(3) * t48 + t33;
t34 = t77 * t46 + t85 * t74;
t47 = -t59 * t74 + t77 * t80;
t32 = qJ(3) * t47 + t34;
t67 = sin(pkin(11));
t70 = cos(pkin(11));
t23 = t67 * t29 + t70 * t32;
t21 = pkin(9) * t55 + t23;
t42 = -t49 * t68 + t71 * t53;
t36 = -pkin(2) * t47 + qJD(3) + t42;
t40 = t47 * t70 - t48 * t67;
t41 = t47 * t67 + t48 * t70;
t25 = -pkin(3) * t40 - pkin(9) * t41 + t36;
t73 = sin(qJ(4));
t83 = cos(qJ(4));
t12 = t83 * t21 + t73 * t25;
t39 = qJD(4) - t40;
t10 = qJ(5) * t39 + t12;
t22 = t29 * t70 - t67 * t32;
t20 = -pkin(3) * t55 - t22;
t37 = t41 * t73 - t83 * t55;
t38 = t41 * t83 + t73 * t55;
t15 = pkin(4) * t37 - qJ(5) * t38 + t20;
t66 = sin(pkin(12));
t69 = cos(pkin(12));
t6 = t69 * t10 + t66 * t15;
t5 = -t10 * t66 + t69 * t15;
t11 = -t73 * t21 + t25 * t83;
t9 = -t39 * pkin(4) + qJD(5) - t11;
t79 = V_base(3) ^ 2;
t76 = cos(qJ(6));
t72 = sin(qJ(6));
t35 = qJD(6) + t37;
t27 = t38 * t69 + t39 * t66;
t26 = -t38 * t66 + t39 * t69;
t17 = t26 * t72 + t27 * t76;
t16 = t26 * t76 - t27 * t72;
t7 = -t26 * pkin(5) + t9;
t4 = pkin(10) * t26 + t6;
t3 = pkin(5) * t37 - pkin(10) * t27 + t5;
t2 = t3 * t72 + t4 * t76;
t1 = t3 * t76 - t4 * t72;
t8 = (-V_base(3) * mrSges(2,1) + t57 * mrSges(2,3) + Ifges(2,4) * t59 + Ifges(2,6) * t65 + Ifges(2,2) * t58 / 0.2e1) * t58 + (t20 * mrSges(5,1) + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t12 * mrSges(5,3) - Ifges(5,4) * t38 + Ifges(6,5) * t27 - Ifges(5,6) * t39 + Ifges(6,6) * t26 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t37) * t37 + (t33 * mrSges(3,1) + t22 * mrSges(4,1) - t34 * mrSges(3,2) - t23 * mrSges(4,2) + Ifges(3,5) * t48 + Ifges(4,5) * t41 + Ifges(3,6) * t47 + Ifges(4,6) * t40 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t55) * t55 + (t56 * mrSges(2,1) - t57 * mrSges(2,2) + Ifges(2,3) * t65 / 0.2e1) * t65 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,2) * t26 / 0.2e1) * t26 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t35 / 0.2e1) * t35 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t39 / 0.2e1) * t39 + (-t36 * mrSges(4,1) + t23 * mrSges(4,3) + Ifges(4,4) * t41 + Ifges(4,2) * t40 / 0.2e1) * t40 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t35 + Ifges(7,1) * t17 / 0.2e1) * t17 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t35 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t27 / 0.2e1) * t27 + (t42 * mrSges(3,2) - t33 * mrSges(3,3) + Ifges(3,1) * t48 / 0.2e1) * t48 + (-t42 * mrSges(3,1) + t34 * mrSges(3,3) + Ifges(3,4) * t48 + Ifges(3,2) * t47 / 0.2e1) * t47 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t20 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t39 + Ifges(5,1) * t38 / 0.2e1) * t38 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t20 ^ 2) / 0.2e1 + (V_base(3) * mrSges(2,2) - t56 * mrSges(2,3) + Ifges(2,5) * t65 + Ifges(2,1) * t59 / 0.2e1) * t59 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(4) * (t22 ^ 2 + t23 ^ 2 + t36 ^ 2) / 0.2e1 + m(3) * (t33 ^ 2 + t34 ^ 2 + t42 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t79) / 0.2e1 + m(2) * (t56 ^ 2 + t57 ^ 2 + t79) / 0.2e1 + (t36 * mrSges(4,2) - t22 * mrSges(4,3) + Ifges(4,1) * t41 / 0.2e1) * t41;
T  = t8;

% Calculate kinetic energy for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR13_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:50
% EndTime: 2019-03-09 14:44:51
% DurationCPUTime: 1.28s
% Computational Cost: add. (3909->160), mult. (5835->227), div. (0->0), fcn. (4850->12), ass. (0->64)
t81 = pkin(2) + pkin(9);
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t54 = t68 * V_base(5) + t72 * V_base(4);
t80 = pkin(8) * t54;
t67 = sin(qJ(2));
t53 = -t68 * V_base(4) + t72 * V_base(5);
t61 = V_base(6) + qJD(1);
t62 = sin(pkin(6));
t63 = cos(pkin(6));
t76 = t53 * t63 + t61 * t62;
t79 = cos(qJ(2));
t42 = t79 * t54 + t76 * t67;
t49 = -t53 * t62 + t61 * t63 + qJD(2);
t59 = V_base(5) * pkin(7) + V_base(1);
t60 = -V_base(4) * pkin(7) + V_base(2);
t51 = t72 * t59 + t68 * t60;
t40 = t76 * pkin(8) + t51;
t50 = -t59 * t68 + t72 * t60;
t43 = pkin(1) * t61 - t63 * t80 + t50;
t46 = -pkin(1) * t53 - t62 * t80 + V_base(3);
t77 = t63 * t79;
t78 = t62 * t79;
t27 = -t67 * t40 + t43 * t77 + t46 * t78;
t74 = qJD(3) - t27;
t19 = t42 * pkin(3) - t81 * t49 + t74;
t41 = -t53 * t77 + t54 * t67 - t61 * t78;
t31 = -t43 * t62 + t63 * t46;
t75 = -qJ(3) * t42 + t31;
t22 = t81 * t41 + t75;
t66 = sin(qJ(4));
t71 = cos(qJ(4));
t12 = t66 * t19 + t71 * t22;
t39 = qJD(4) + t42;
t10 = pkin(10) * t39 + t12;
t28 = t79 * t40 + (t43 * t63 + t46 * t62) * t67;
t26 = -t49 * qJ(3) - t28;
t23 = -pkin(3) * t41 - t26;
t34 = t41 * t71 - t66 * t49;
t35 = t41 * t66 + t49 * t71;
t15 = -pkin(4) * t34 - pkin(10) * t35 + t23;
t65 = sin(qJ(5));
t70 = cos(qJ(5));
t6 = t70 * t10 + t65 * t15;
t5 = -t10 * t65 + t70 * t15;
t11 = t19 * t71 - t66 * t22;
t33 = qJD(5) - t34;
t9 = -pkin(4) * t39 - t11;
t73 = V_base(3) ^ 2;
t69 = cos(qJ(6));
t64 = sin(qJ(6));
t32 = qJD(6) + t33;
t30 = t35 * t70 + t39 * t65;
t29 = -t35 * t65 + t39 * t70;
t25 = -t49 * pkin(2) + t74;
t24 = pkin(2) * t41 + t75;
t17 = t29 * t64 + t30 * t69;
t16 = t29 * t69 - t30 * t64;
t7 = -pkin(5) * t29 + t9;
t4 = pkin(11) * t29 + t6;
t3 = pkin(5) * t33 - pkin(11) * t30 + t5;
t2 = t3 * t64 + t4 * t69;
t1 = t3 * t69 - t4 * t64;
t8 = m(2) * (t50 ^ 2 + t51 ^ 2 + t73) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t73) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(3) * (t27 ^ 2 + t28 ^ 2 + t31 ^ 2) / 0.2e1 + (V_base(3) * mrSges(2,2) - t50 * mrSges(2,3) + Ifges(2,5) * t61 + Ifges(2,1) * t54 / 0.2e1) * t54 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t32 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t31 * mrSges(3,1) + t26 * mrSges(4,1) - t24 * mrSges(4,2) - t28 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t41 + (Ifges(4,5) - Ifges(3,6)) * t49 + (-Ifges(3,4) - Ifges(4,6)) * t42) * t41 + (t25 * mrSges(4,1) + t31 * mrSges(3,2) - t27 * mrSges(3,3) - t24 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t42 + (-Ifges(4,4) + Ifges(3,5)) * t49) * t42 + (-V_base(3) * mrSges(2,1) + t51 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,6) * t61 + Ifges(2,2) * t53 / 0.2e1) * t53 + m(5) * (t11 ^ 2 + t12 ^ 2 + t23 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (t23 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t39 + Ifges(5,1) * t35 / 0.2e1) * t35 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t32 / 0.2e1) * t32 + (-t23 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t35 + Ifges(5,6) * t39 + Ifges(5,2) * t34 / 0.2e1) * t34 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t33 / 0.2e1) * t33 + (t27 * mrSges(3,1) - t28 * mrSges(3,2) + t25 * mrSges(4,2) - t26 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t49) * t49 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t30 + Ifges(6,6) * t33 + Ifges(6,2) * t29 / 0.2e1) * t29 + (t50 * mrSges(2,1) - t51 * mrSges(2,2) + Ifges(2,3) * t61 / 0.2e1) * t61 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t32 + Ifges(7,2) * t16 / 0.2e1) * t16 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t33 + Ifges(6,1) * t30 / 0.2e1) * t30 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t39 / 0.2e1) * t39 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5);
T  = t8;

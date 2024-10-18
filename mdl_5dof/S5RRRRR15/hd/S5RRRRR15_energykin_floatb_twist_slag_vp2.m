% Calculate kinetic energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR15_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:52
% EndTime: 2024-09-27 22:23:53
% DurationCPUTime: 0.44s
% Computational Cost: add. (4807->142), mult. (7724->215), div. (0->0), fcn. (6602->14), ass. (0->63)
t55 = V_base(5) * pkin(7) + V_base(1);
t56 = -V_base(4) * pkin(7) + V_base(2);
t66 = sin(qJ(1));
t71 = cos(qJ(1));
t47 = -t55 * t66 + t71 * t56;
t57 = V_base(6) + qJD(1);
t61 = cos(pkin(5));
t50 = t66 * V_base(5) + t71 * V_base(4);
t79 = pkin(8) * t50;
t39 = pkin(1) * t57 - t61 * t79 + t47;
t49 = -t66 * V_base(4) + t71 * V_base(5);
t59 = sin(pkin(5));
t43 = -pkin(1) * t49 - t59 * t79 + V_base(3);
t80 = t39 * t61 + t43 * t59;
t65 = sin(qJ(2));
t70 = cos(qJ(2));
t73 = t49 * t61 + t57 * t59;
t37 = -t50 * t65 + t73 * t70;
t38 = t50 * t70 + t73 * t65;
t64 = sin(qJ(3));
t69 = cos(qJ(3));
t30 = t37 * t69 - t38 * t64;
t31 = t37 * t64 + t38 * t69;
t63 = sin(qJ(4));
t68 = cos(qJ(4));
t24 = t30 * t63 + t31 * t68;
t78 = pkin(11) * t24;
t48 = t71 * t55 + t66 * t56;
t36 = t73 * pkin(8) + t48;
t27 = -t36 * t65 + t80 * t70;
t46 = -t49 * t59 + t61 * t57 + qJD(2);
t22 = pkin(2) * t46 - pkin(9) * t38 + t27;
t28 = t70 * t36 + t80 * t65;
t26 = pkin(9) * t37 + t28;
t16 = t69 * t22 - t26 * t64;
t45 = qJD(3) + t46;
t11 = pkin(3) * t45 - pkin(10) * t31 + t16;
t17 = t64 * t22 + t69 * t26;
t13 = pkin(10) * t30 + t17;
t7 = t63 * t11 + t68 * t13;
t6 = t68 * t11 - t13 * t63;
t32 = -t39 * t59 + t61 * t43;
t44 = qJD(4) + t45;
t60 = cos(pkin(6));
t5 = pkin(4) * t44 - t60 * t78 + t6;
t58 = sin(pkin(6));
t29 = -pkin(2) * t37 + t32;
t19 = -pkin(3) * t30 + t29;
t23 = t30 * t68 - t31 * t63;
t8 = -pkin(4) * t23 - t58 * t78 + t19;
t75 = t5 * t60 + t58 * t8;
t74 = t23 * t60 + t44 * t58;
t72 = V_base(3) ^ 2;
t67 = cos(qJ(5));
t62 = sin(qJ(5));
t18 = -t23 * t58 + t44 * t60 + qJD(5);
t15 = t24 * t67 + t74 * t62;
t14 = -t24 * t62 + t74 * t67;
t4 = t74 * pkin(11) + t7;
t3 = -t5 * t58 + t60 * t8;
t2 = t4 * t67 + t75 * t62;
t1 = -t4 * t62 + t75 * t67;
t9 = m(4) * (t16 ^ 2 + t17 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t27 ^ 2 + t28 ^ 2 + t32 ^ 2) / 0.2e1 + m(2) * (t47 ^ 2 + t48 ^ 2 + t72) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t72) / 0.2e1 + m(5) * (t19 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t47 * mrSges(2,1) - t48 * mrSges(2,2) + Ifges(2,3) * t57 / 0.2e1) * t57 + (t27 * mrSges(3,1) - t28 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t45 / 0.2e1) * t45 + (t6 * mrSges(5,1) - t7 * mrSges(5,2) + Ifges(5,3) * t44 / 0.2e1) * t44 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t18 / 0.2e1) * t18 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t47 * mrSges(2,3) + Ifges(2,5) * t57 + Ifges(2,1) * t50 / 0.2e1) * t50 + (t32 * mrSges(3,2) - t27 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t29 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t45 + Ifges(4,1) * t31 / 0.2e1) * t31 + (t19 * mrSges(5,2) - t6 * mrSges(5,3) + Ifges(5,5) * t44 + Ifges(5,1) * t24 / 0.2e1) * t24 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t18 + Ifges(6,1) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t48 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,6) * t57 + Ifges(2,2) * t49 / 0.2e1) * t49 + (-t32 * mrSges(3,1) + t28 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t46 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t29 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t31 + Ifges(4,6) * t45 + Ifges(4,2) * t30 / 0.2e1) * t30 + (-t19 * mrSges(5,1) + t7 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,6) * t44 + Ifges(5,2) * t23 / 0.2e1) * t23 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t18 + Ifges(6,2) * t14 / 0.2e1) * t14;
T = t9;

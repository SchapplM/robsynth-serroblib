% Calculate kinetic energy for
% S5PRRRR12
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:08
% EndTime: 2024-09-28 18:07:09
% DurationCPUTime: 0.55s
% Computational Cost: add. (4205->142), mult. (7724->214), div. (0->0), fcn. (6602->14), ass. (0->62)
t54 = V_base(5) * qJ(1) + V_base(1);
t55 = -V_base(4) * qJ(1) + V_base(2);
t58 = sin(pkin(11));
t61 = cos(pkin(11));
t47 = -t54 * t58 + t61 * t55;
t63 = cos(pkin(5));
t50 = t58 * V_base(5) + t61 * V_base(4);
t78 = pkin(7) * t50;
t39 = V_base(6) * pkin(1) - t63 * t78 + t47;
t49 = -t58 * V_base(4) + t61 * V_base(5);
t57 = V_base(3) + qJD(1);
t60 = sin(pkin(5));
t43 = -pkin(1) * t49 - t60 * t78 + t57;
t79 = t39 * t63 + t43 * t60;
t67 = sin(qJ(2));
t71 = cos(qJ(2));
t72 = t49 * t63 + t60 * V_base(6);
t37 = -t67 * t50 + t72 * t71;
t38 = t50 * t71 + t72 * t67;
t66 = sin(qJ(3));
t70 = cos(qJ(3));
t30 = t37 * t70 - t38 * t66;
t31 = t37 * t66 + t38 * t70;
t65 = sin(qJ(4));
t69 = cos(qJ(4));
t24 = t30 * t65 + t31 * t69;
t77 = pkin(10) * t24;
t48 = t61 * t54 + t58 * t55;
t36 = t72 * pkin(7) + t48;
t27 = -t36 * t67 + t79 * t71;
t46 = -t49 * t60 + t63 * V_base(6) + qJD(2);
t22 = pkin(2) * t46 - pkin(8) * t38 + t27;
t28 = t71 * t36 + t79 * t67;
t26 = pkin(8) * t37 + t28;
t16 = t70 * t22 - t26 * t66;
t45 = qJD(3) + t46;
t11 = pkin(3) * t45 - pkin(9) * t31 + t16;
t17 = t66 * t22 + t70 * t26;
t13 = pkin(9) * t30 + t17;
t7 = t65 * t11 + t69 * t13;
t6 = t69 * t11 - t13 * t65;
t32 = -t39 * t60 + t63 * t43;
t44 = qJD(4) + t45;
t62 = cos(pkin(6));
t5 = pkin(4) * t44 - t62 * t77 + t6;
t59 = sin(pkin(6));
t29 = -pkin(2) * t37 + t32;
t19 = -pkin(3) * t30 + t29;
t23 = t30 * t69 - t31 * t65;
t8 = -pkin(4) * t23 - t59 * t77 + t19;
t74 = t5 * t62 + t59 * t8;
t73 = t23 * t62 + t44 * t59;
t68 = cos(qJ(5));
t64 = sin(qJ(5));
t18 = -t23 * t59 + t44 * t62 + qJD(5);
t15 = t24 * t68 + t73 * t64;
t14 = -t24 * t64 + t73 * t68;
t4 = t73 * pkin(10) + t7;
t3 = -t5 * t59 + t62 * t8;
t2 = t4 * t68 + t74 * t64;
t1 = -t4 * t64 + t74 * t68;
t9 = m(2) * (t47 ^ 2 + t48 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t29 ^ 2) / 0.2e1 + m(3) * (t27 ^ 2 + t28 ^ 2 + t32 ^ 2) / 0.2e1 + m(5) * (t19 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t57 * mrSges(2,2) - t47 * mrSges(2,3) + Ifges(2,1) * t50 / 0.2e1) * t50 + (t27 * mrSges(3,1) - t28 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + Ifges(4,3) * t45 / 0.2e1) * t45 + (t6 * mrSges(5,1) - t7 * mrSges(5,2) + Ifges(5,3) * t44 / 0.2e1) * t44 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t18 / 0.2e1) * t18 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t57 * mrSges(2,1) + t48 * mrSges(2,3) + Ifges(2,4) * t50 + Ifges(2,2) * t49 / 0.2e1) * t49 + (t32 * mrSges(3,2) - t27 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t29 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,5) * t45 + Ifges(4,1) * t31 / 0.2e1) * t31 + (t19 * mrSges(5,2) - t6 * mrSges(5,3) + Ifges(5,5) * t44 + Ifges(5,1) * t24 / 0.2e1) * t24 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t18 + Ifges(6,1) * t15 / 0.2e1) * t15 + (-t32 * mrSges(3,1) + t28 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t46 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t29 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t31 + Ifges(4,6) * t45 + Ifges(4,2) * t30 / 0.2e1) * t30 + (-t19 * mrSges(5,1) + t7 * mrSges(5,3) + Ifges(5,4) * t24 + Ifges(5,6) * t44 + Ifges(5,2) * t23 / 0.2e1) * t23 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t18 + Ifges(6,2) * t14 / 0.2e1) * t14 + (V_base(2) * mrSges(1,1) + t47 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t48 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t50 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t49 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t9;

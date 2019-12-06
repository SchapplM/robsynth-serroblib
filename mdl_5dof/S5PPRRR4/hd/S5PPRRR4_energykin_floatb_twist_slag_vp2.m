% Calculate kinetic energy for
% S5PPRRR4
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:31
% EndTime: 2019-12-05 15:18:32
% DurationCPUTime: 1.25s
% Computational Cost: add. (4239->142), mult. (7958->214), div. (0->0), fcn. (6794->14), ass. (0->62)
t55 = V_base(5) * qJ(1) + V_base(1);
t56 = -V_base(4) * qJ(1) + V_base(2);
t59 = sin(pkin(10));
t63 = cos(pkin(10));
t48 = -t55 * t59 + t63 * t56;
t65 = cos(pkin(5));
t51 = t59 * V_base(5) + t63 * V_base(4);
t75 = qJ(2) * t51;
t42 = V_base(6) * pkin(1) - t65 * t75 + t48;
t50 = -t59 * V_base(4) + t63 * V_base(5);
t57 = V_base(3) + qJD(1);
t61 = sin(pkin(5));
t46 = -pkin(1) * t50 - t61 * t75 + t57;
t83 = t42 * t65 + t46 * t61;
t49 = t63 * t55 + t59 * t56;
t72 = t50 * t65 + t61 * V_base(6);
t39 = t72 * qJ(2) + t49;
t58 = sin(pkin(11));
t62 = cos(pkin(11));
t29 = -t39 * t58 + t83 * t62;
t47 = -t50 * t61 + t65 * V_base(6);
t64 = cos(pkin(6));
t41 = t51 * t62 + t72 * t58;
t81 = pkin(7) * t41;
t22 = pkin(2) * t47 - t64 * t81 + t29;
t34 = -t42 * t61 + t65 * t46 + qJD(2);
t40 = -t51 * t58 + t72 * t62;
t60 = sin(pkin(6));
t28 = -pkin(2) * t40 - t60 * t81 + t34;
t82 = t22 * t64 + t28 * t60;
t30 = t62 * t39 + t83 * t58;
t73 = t40 * t64 + t47 * t60;
t21 = t73 * pkin(7) + t30;
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t13 = -t68 * t21 + t82 * t71;
t32 = -t68 * t41 + t73 * t71;
t14 = t71 * t21 + t82 * t68;
t35 = -t40 * t60 + t47 * t64 + qJD(3);
t10 = pkin(8) * t35 + t14;
t15 = -t22 * t60 + t64 * t28;
t33 = t41 * t71 + t73 * t68;
t12 = -pkin(3) * t32 - pkin(8) * t33 + t15;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t6 = t70 * t10 + t67 * t12;
t5 = -t10 * t67 + t12 * t70;
t24 = -t33 * t67 + t35 * t70;
t9 = -t35 * pkin(3) - t13;
t69 = cos(qJ(5));
t66 = sin(qJ(5));
t31 = qJD(4) - t32;
t25 = t33 * t70 + t35 * t67;
t23 = qJD(5) - t24;
t17 = t25 * t69 + t31 * t66;
t16 = -t25 * t66 + t31 * t69;
t7 = -t24 * pkin(4) - t25 * pkin(9) + t9;
t4 = pkin(9) * t31 + t6;
t3 = -pkin(4) * t31 - t5;
t2 = t4 * t69 + t66 * t7;
t1 = -t4 * t66 + t69 * t7;
t8 = m(4) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t34 ^ 2) / 0.2e1 + m(2) * (t48 ^ 2 + t49 ^ 2 + t57 ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t57 * mrSges(2,2) - t48 * mrSges(2,3) + Ifges(2,1) * t51 / 0.2e1) * t51 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t47 / 0.2e1) * t47 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t35 / 0.2e1) * t35 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t31 / 0.2e1) * t31 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t57 * mrSges(2,1) + t49 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,2) * t50 / 0.2e1) * t50 + (t34 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t47 + Ifges(3,1) * t41 / 0.2e1) * t41 + (t15 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t35 + Ifges(4,1) * t33 / 0.2e1) * t33 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t31 + Ifges(5,1) * t25 / 0.2e1) * t25 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t23 + Ifges(6,1) * t17 / 0.2e1) * t17 + (-t34 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t41 + Ifges(3,6) * t47 + Ifges(3,2) * t40 / 0.2e1) * t40 + (-t15 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t33 + Ifges(4,6) * t35 + Ifges(4,2) * t32 / 0.2e1) * t32 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,6) * t31 + Ifges(5,2) * t24 / 0.2e1) * t24 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t23 + Ifges(6,2) * t16 / 0.2e1) * t16 + (V_base(2) * mrSges(1,1) + t48 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t49 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t51 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t50 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t8;

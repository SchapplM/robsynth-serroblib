% Calculate kinetic energy for
% S5RRRRR12
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:29
% EndTime: 2019-12-31 22:46:30
% DurationCPUTime: 1.32s
% Computational Cost: add. (5065->142), mult. (7958->215), div. (0->0), fcn. (6794->14), ass. (0->63)
t55 = V_base(5) * pkin(7) + V_base(1);
t56 = -V_base(4) * pkin(7) + V_base(2);
t66 = sin(qJ(1));
t71 = cos(qJ(1));
t48 = -t55 * t66 + t71 * t56;
t57 = V_base(6) + qJD(1);
t61 = cos(pkin(5));
t51 = t66 * V_base(5) + t71 * V_base(4);
t82 = pkin(8) * t51;
t42 = pkin(1) * t57 - t61 * t82 + t48;
t50 = -t66 * V_base(4) + t71 * V_base(5);
t59 = sin(pkin(5));
t46 = -pkin(1) * t50 - t59 * t82 + V_base(3);
t84 = t42 * t61 + t46 * t59;
t49 = t71 * t55 + t66 * t56;
t73 = t50 * t61 + t57 * t59;
t39 = t73 * pkin(8) + t49;
t65 = sin(qJ(2));
t70 = cos(qJ(2));
t29 = -t39 * t65 + t84 * t70;
t47 = -t50 * t59 + t57 * t61 + qJD(2);
t60 = cos(pkin(6));
t41 = t51 * t70 + t73 * t65;
t81 = pkin(9) * t41;
t25 = pkin(2) * t47 - t60 * t81 + t29;
t34 = -t42 * t59 + t61 * t46;
t40 = -t51 * t65 + t73 * t70;
t58 = sin(pkin(6));
t28 = -pkin(2) * t40 - t58 * t81 + t34;
t83 = t25 * t60 + t28 * t58;
t30 = t70 * t39 + t84 * t65;
t74 = t40 * t60 + t47 * t58;
t24 = t74 * pkin(9) + t30;
t64 = sin(qJ(3));
t69 = cos(qJ(3));
t13 = -t64 * t24 + t83 * t69;
t32 = -t41 * t64 + t74 * t69;
t14 = t69 * t24 + t83 * t64;
t35 = -t40 * t58 + t47 * t60 + qJD(3);
t10 = pkin(10) * t35 + t14;
t15 = -t25 * t58 + t60 * t28;
t33 = t41 * t69 + t74 * t64;
t12 = -pkin(3) * t32 - pkin(10) * t33 + t15;
t63 = sin(qJ(4));
t68 = cos(qJ(4));
t6 = t68 * t10 + t63 * t12;
t5 = -t10 * t63 + t12 * t68;
t22 = -t33 * t63 + t35 * t68;
t9 = -pkin(3) * t35 - t13;
t72 = V_base(3) ^ 2;
t67 = cos(qJ(5));
t62 = sin(qJ(5));
t31 = qJD(4) - t32;
t23 = t33 * t68 + t35 * t63;
t21 = qJD(5) - t22;
t17 = t23 * t67 + t31 * t62;
t16 = -t23 * t62 + t31 * t67;
t7 = -pkin(4) * t22 - pkin(11) * t23 + t9;
t4 = pkin(11) * t31 + t6;
t3 = -pkin(4) * t31 - t5;
t2 = t4 * t67 + t62 * t7;
t1 = -t4 * t62 + t67 * t7;
t8 = m(2) * (t48 ^ 2 + t49 ^ 2 + t72) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t72) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t34 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t48 * mrSges(2,1) - t49 * mrSges(2,2) + Ifges(2,3) * t57 / 0.2e1) * t57 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t47 / 0.2e1) * t47 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t35 / 0.2e1) * t35 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t31 / 0.2e1) * t31 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t21 / 0.2e1) * t21 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t48 * mrSges(2,3) + Ifges(2,5) * t57 + Ifges(2,1) * t51 / 0.2e1) * t51 + (t34 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t47 + Ifges(3,1) * t41 / 0.2e1) * t41 + (t15 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t35 + Ifges(4,1) * t33 / 0.2e1) * t33 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t31 + Ifges(5,1) * t23 / 0.2e1) * t23 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t21 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t49 * mrSges(2,3) + Ifges(2,4) * t51 + Ifges(2,6) * t57 + Ifges(2,2) * t50 / 0.2e1) * t50 + (-t34 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t41 + Ifges(3,6) * t47 + Ifges(3,2) * t40 / 0.2e1) * t40 + (-t15 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t33 + Ifges(4,6) * t35 + Ifges(4,2) * t32 / 0.2e1) * t32 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t23 + Ifges(5,6) * t31 + Ifges(5,2) * t22 / 0.2e1) * t22 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t21 + Ifges(6,2) * t16 / 0.2e1) * t16;
T = t8;

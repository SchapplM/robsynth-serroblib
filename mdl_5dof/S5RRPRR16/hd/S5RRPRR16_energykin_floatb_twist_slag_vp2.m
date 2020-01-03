% Calculate kinetic energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR16_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:37
% EndTime: 2019-12-31 20:44:38
% DurationCPUTime: 0.91s
% Computational Cost: add. (1985->136), mult. (2976->192), div. (0->0), fcn. (2376->10), ass. (0->54)
t66 = pkin(2) + pkin(8);
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t41 = t54 * V_base(5) + t57 * V_base(4);
t65 = pkin(7) * t41;
t40 = -t54 * V_base(4) + t57 * V_base(5);
t48 = V_base(6) + qJD(1);
t53 = sin(qJ(2));
t50 = cos(pkin(5));
t64 = cos(qJ(2));
t62 = t50 * t64;
t49 = sin(pkin(5));
t63 = t49 * t64;
t29 = -t40 * t62 + t41 * t53 - t48 * t63;
t46 = V_base(5) * pkin(6) + V_base(1);
t47 = -V_base(4) * pkin(6) + V_base(2);
t37 = -t46 * t54 + t57 * t47;
t31 = pkin(1) * t48 - t50 * t65 + t37;
t34 = -pkin(1) * t40 - t49 * t65 + V_base(3);
t20 = -t31 * t49 + t50 * t34;
t61 = t40 * t50 + t48 * t49;
t30 = t64 * t41 + t61 * t53;
t60 = -qJ(3) * t30 + t20;
t11 = t66 * t29 + t60;
t52 = sin(qJ(4));
t56 = cos(qJ(4));
t36 = -t40 * t49 + t48 * t50 + qJD(2);
t38 = t57 * t46 + t54 * t47;
t28 = t61 * pkin(7) + t38;
t16 = -t53 * t28 + t31 * t62 + t34 * t63;
t59 = qJD(3) - t16;
t9 = t30 * pkin(3) - t66 * t36 + t59;
t6 = t56 * t11 + t52 * t9;
t17 = t64 * t28 + (t31 * t50 + t34 * t49) * t53;
t15 = -t36 * qJ(3) - t17;
t5 = -t11 * t52 + t56 * t9;
t22 = t29 * t56 - t36 * t52;
t12 = -pkin(3) * t29 - t15;
t58 = V_base(3) ^ 2;
t55 = cos(qJ(5));
t51 = sin(qJ(5));
t27 = qJD(4) + t30;
t23 = t29 * t52 + t36 * t56;
t21 = qJD(5) - t22;
t19 = t23 * t55 + t27 * t51;
t18 = -t23 * t51 + t27 * t55;
t14 = -t36 * pkin(2) + t59;
t13 = pkin(2) * t29 + t60;
t7 = -pkin(4) * t22 - pkin(9) * t23 + t12;
t4 = pkin(9) * t27 + t6;
t3 = -pkin(4) * t27 - t5;
t2 = t4 * t55 + t51 * t7;
t1 = -t4 * t51 + t55 * t7;
t8 = m(2) * (t37 ^ 2 + t38 ^ 2 + t58) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t58) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(3) * (t16 ^ 2 + t17 ^ 2 + t20 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t37 * mrSges(2,1) - t38 * mrSges(2,2) + Ifges(2,3) * t48 / 0.2e1) * t48 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t27 / 0.2e1) * t27 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t21 / 0.2e1) * t21 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t37 * mrSges(2,3) + Ifges(2,5) * t48 + Ifges(2,1) * t41 / 0.2e1) * t41 + (t12 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t27 + Ifges(5,1) * t23 / 0.2e1) * t23 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t21 + Ifges(6,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t38 * mrSges(2,3) + Ifges(2,4) * t41 + Ifges(2,6) * t48 + Ifges(2,2) * t40 / 0.2e1) * t40 + (-t12 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t23 + Ifges(5,6) * t27 + Ifges(5,2) * t22 / 0.2e1) * t22 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t19 + Ifges(6,6) * t21 + Ifges(6,2) * t18 / 0.2e1) * t18 + (t16 * mrSges(3,1) - t17 * mrSges(3,2) + t14 * mrSges(4,2) - t15 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t36) * t36 + (t14 * mrSges(4,1) + t20 * mrSges(3,2) - t16 * mrSges(3,3) - t13 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t30 + (-Ifges(4,4) + Ifges(3,5)) * t36) * t30 + (t20 * mrSges(3,1) + t15 * mrSges(4,1) - t13 * mrSges(4,2) - t17 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t29 + (Ifges(4,5) - Ifges(3,6)) * t36 + (-Ifges(3,4) - Ifges(4,6)) * t30) * t29;
T = t8;

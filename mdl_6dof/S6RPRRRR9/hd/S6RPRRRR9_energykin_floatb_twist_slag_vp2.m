% Calculate kinetic energy for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:01
% EndTime: 2019-03-09 07:23:02
% DurationCPUTime: 1.08s
% Computational Cost: add. (2539->153), mult. (3191->213), div. (0->0), fcn. (2396->10), ass. (0->57)
t71 = pkin(1) + pkin(7);
t62 = sin(qJ(1));
t70 = cos(qJ(1));
t48 = t62 * V_base(5) + t70 * V_base(4);
t57 = V_base(6) + qJD(1);
t53 = V_base(5) * pkin(6) + V_base(1);
t54 = -V_base(4) * pkin(6) + V_base(2);
t44 = -t62 * t53 + t54 * t70;
t68 = qJD(2) - t44;
t29 = t48 * pkin(2) - t57 * t71 + t68;
t47 = t62 * V_base(4) - t70 * V_base(5);
t69 = -qJ(2) * t48 + V_base(3);
t34 = t47 * t71 + t69;
t61 = sin(qJ(3));
t66 = cos(qJ(3));
t22 = t61 * t29 + t66 * t34;
t46 = qJD(3) + t48;
t20 = pkin(8) * t46 + t22;
t45 = t70 * t53 + t62 * t54;
t39 = -t57 * qJ(2) - t45;
t35 = -pkin(2) * t47 - t39;
t42 = t47 * t66 - t61 * t57;
t43 = t47 * t61 + t57 * t66;
t27 = -pkin(3) * t42 - pkin(8) * t43 + t35;
t60 = sin(qJ(4));
t65 = cos(qJ(4));
t14 = t65 * t20 + t60 * t27;
t32 = -t43 * t60 + t46 * t65;
t11 = pkin(9) * t32 + t14;
t59 = sin(qJ(5));
t64 = cos(qJ(5));
t13 = -t20 * t60 + t65 * t27;
t33 = t43 * t65 + t46 * t60;
t41 = qJD(4) - t42;
t9 = pkin(4) * t41 - pkin(9) * t33 + t13;
t6 = t64 * t11 + t59 * t9;
t5 = -t11 * t59 + t64 * t9;
t21 = t29 * t66 - t61 * t34;
t19 = -pkin(3) * t46 - t21;
t40 = qJD(5) + t41;
t17 = -pkin(4) * t32 + t19;
t67 = V_base(3) ^ 2;
t63 = cos(qJ(6));
t58 = sin(qJ(6));
t38 = -t57 * pkin(1) + t68;
t37 = qJD(6) + t40;
t36 = pkin(1) * t47 + t69;
t24 = t32 * t59 + t33 * t64;
t23 = t32 * t64 - t33 * t59;
t16 = t23 * t58 + t24 * t63;
t15 = t23 * t63 - t24 * t58;
t12 = -pkin(5) * t23 + t17;
t4 = pkin(10) * t23 + t6;
t3 = pkin(5) * t40 - pkin(10) * t24 + t5;
t2 = t3 * t58 + t4 * t63;
t1 = t3 * t63 - t4 * t58;
t7 = (-t17 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t24 + Ifges(6,6) * t40 + Ifges(6,2) * t23 / 0.2e1) * t23 + (t13 * mrSges(5,1) - t14 * mrSges(5,2) + Ifges(5,3) * t41 / 0.2e1) * t41 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t12 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t37 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t35 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t46 + Ifges(4,1) * t43 / 0.2e1) * t43 + (-t35 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t46 + Ifges(4,2) * t42 / 0.2e1) * t42 + (t12 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t37 + Ifges(7,1) * t16 / 0.2e1) * t16 + m(2) * (t44 ^ 2 + t45 ^ 2 + t67) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t67) / 0.2e1 + m(3) * (t36 ^ 2 + t38 ^ 2 + t39 ^ 2) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t35 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) / 0.2e1 + m(6) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t17 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t40 + Ifges(6,1) * t24 / 0.2e1) * t24 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t40 / 0.2e1) * t40 + (V_base(3) * mrSges(2,1) + t39 * mrSges(3,1) - t36 * mrSges(3,2) - t45 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t47 + (Ifges(3,5) - Ifges(2,6)) * t57 + (-Ifges(2,4) - Ifges(3,6)) * t48) * t47 + (t38 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) - t36 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t48 + (-Ifges(3,4) + Ifges(2,5)) * t57) * t48 + (t19 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,5) * t41 + Ifges(5,1) * t33 / 0.2e1) * t33 + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t46 / 0.2e1) * t46 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t37 / 0.2e1) * t37 + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + t38 * mrSges(3,2) - t39 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t57) * t57 + (-t19 * mrSges(5,1) + t14 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,6) * t41 + Ifges(5,2) * t32 / 0.2e1) * t32;
T  = t7;

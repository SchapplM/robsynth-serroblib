% Calculate kinetic energy for
% S6RPRRRR8
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:35
% EndTime: 2019-03-09 07:19:36
% DurationCPUTime: 1.08s
% Computational Cost: add. (2447->153), mult. (3055->213), div. (0->0), fcn. (2284->10), ass. (0->57)
t71 = pkin(1) + pkin(7);
t62 = sin(qJ(1));
t67 = cos(qJ(1));
t49 = t62 * V_base(5) + t67 * V_base(4);
t57 = V_base(6) + qJD(1);
t53 = V_base(5) * pkin(6) + V_base(1);
t54 = -V_base(4) * pkin(6) + V_base(2);
t44 = -t62 * t53 + t54 * t67;
t69 = qJD(2) - t44;
t34 = pkin(2) * t49 - t57 * t71 + t69;
t48 = t62 * V_base(4) - t67 * V_base(5);
t70 = -qJ(2) * t49 + V_base(3);
t36 = t48 * t71 + t70;
t61 = sin(qJ(3));
t66 = cos(qJ(3));
t23 = t66 * t34 - t36 * t61;
t43 = t48 * t61 + t57 * t66;
t47 = qJD(3) + t49;
t19 = pkin(3) * t47 - pkin(8) * t43 + t23;
t24 = t61 * t34 + t66 * t36;
t42 = t48 * t66 - t57 * t61;
t22 = pkin(8) * t42 + t24;
t60 = sin(qJ(4));
t65 = cos(qJ(4));
t12 = t60 * t19 + t65 * t22;
t46 = qJD(4) + t47;
t10 = pkin(9) * t46 + t12;
t45 = t67 * t53 + t62 * t54;
t41 = -t57 * qJ(2) - t45;
t37 = -pkin(2) * t48 - t41;
t27 = -pkin(3) * t42 + t37;
t30 = t42 * t65 - t60 * t43;
t31 = t42 * t60 + t43 * t65;
t15 = -pkin(4) * t30 - pkin(9) * t31 + t27;
t59 = sin(qJ(5));
t64 = cos(qJ(5));
t6 = t64 * t10 + t59 * t15;
t5 = -t10 * t59 + t64 * t15;
t11 = t19 * t65 - t60 * t22;
t29 = qJD(5) - t30;
t9 = -pkin(4) * t46 - t11;
t68 = V_base(3) ^ 2;
t63 = cos(qJ(6));
t58 = sin(qJ(6));
t39 = -pkin(1) * t57 + t69;
t38 = pkin(1) * t48 + t70;
t28 = qJD(6) + t29;
t26 = t31 * t64 + t46 * t59;
t25 = -t31 * t59 + t46 * t64;
t17 = t25 * t58 + t26 * t63;
t16 = t25 * t63 - t26 * t58;
t7 = -pkin(5) * t25 + t9;
t4 = pkin(10) * t25 + t6;
t3 = pkin(5) * t29 - pkin(10) * t26 + t5;
t2 = t3 * t58 + t4 * t63;
t1 = t3 * t63 - t4 * t58;
t8 = (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t26 / 0.2e1) * t26 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t46 / 0.2e1) * t46 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t29 + Ifges(6,2) * t25 / 0.2e1) * t25 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t28 + Ifges(7,2) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(2,1) + t41 * mrSges(3,1) - t38 * mrSges(3,2) - t45 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t48 + (Ifges(3,5) - Ifges(2,6)) * t57 + (-Ifges(2,4) - Ifges(3,6)) * t49) * t48 + (t39 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) - t38 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(2,1) / 0.2e1) * t49 + (-Ifges(3,4) + Ifges(2,5)) * t57) * t49 + (t27 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t46 + Ifges(5,1) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t37 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,5) * t47 + Ifges(4,1) * t43 / 0.2e1) * t43 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t28 / 0.2e1) * t28 + m(2) * (t44 ^ 2 + t45 ^ 2 + t68) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t68) / 0.2e1 + m(3) * (t38 ^ 2 + t39 ^ 2 + t41 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t28 + Ifges(7,1) * t17 / 0.2e1) * t17 + (t23 * mrSges(4,1) - t24 * mrSges(4,2) + Ifges(4,3) * t47 / 0.2e1) * t47 + m(4) * (t23 ^ 2 + t24 ^ 2 + t37 ^ 2) / 0.2e1 + (-t37 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t43 + Ifges(4,6) * t47 + Ifges(4,2) * t42 / 0.2e1) * t42 + m(5) * (t11 ^ 2 + t12 ^ 2 + t27 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + t39 * mrSges(3,2) - t41 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t57) * t57 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (-t27 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t31 + Ifges(5,6) * t46 + Ifges(5,2) * t30 / 0.2e1) * t30;
T  = t8;

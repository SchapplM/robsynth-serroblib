% Calculate kinetic energy for
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:33:59
% EndTime: 2019-03-09 09:34:00
% DurationCPUTime: 1.00s
% Computational Cost: add. (2625->153), mult. (3327->211), div. (0->0), fcn. (2568->10), ass. (0->56)
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t48 = t61 * V_base(5) + t65 * V_base(4);
t55 = V_base(6) + qJD(1);
t60 = sin(qJ(2));
t64 = cos(qJ(2));
t42 = t48 * t64 + t55 * t60;
t47 = -t61 * V_base(4) + t65 * V_base(5);
t46 = qJD(2) - t47;
t34 = -pkin(1) * t47 - pkin(7) * t48 + V_base(3);
t53 = V_base(5) * pkin(6) + V_base(1);
t54 = -V_base(4) * pkin(6) + V_base(2);
t44 = t65 * t53 + t61 * t54;
t39 = pkin(7) * t55 + t44;
t29 = t34 * t64 - t60 * t39;
t68 = qJD(3) - t29;
t69 = pkin(2) + qJ(4);
t20 = pkin(3) * t42 - t69 * t46 + t68;
t41 = t48 * t60 - t64 * t55;
t43 = -t61 * t53 + t54 * t65;
t38 = -pkin(1) * t55 - t43;
t67 = -qJ(3) * t42 + t38;
t23 = t69 * t41 + t67;
t56 = sin(pkin(10));
t57 = cos(pkin(10));
t14 = t56 * t20 + t57 * t23;
t31 = t41 * t57 - t46 * t56;
t11 = pkin(8) * t31 + t14;
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t13 = t57 * t20 - t23 * t56;
t32 = t41 * t56 + t46 * t57;
t9 = pkin(4) * t42 - pkin(8) * t32 + t13;
t6 = t63 * t11 + t59 * t9;
t30 = t60 * t34 + t64 * t39;
t27 = -t46 * qJ(3) - t30;
t5 = -t11 * t59 + t63 * t9;
t40 = qJD(5) + t42;
t22 = -pkin(3) * t41 + qJD(4) - t27;
t17 = -pkin(4) * t31 + t22;
t66 = V_base(3) ^ 2;
t62 = cos(qJ(6));
t58 = sin(qJ(6));
t37 = qJD(6) + t40;
t28 = pkin(2) * t41 + t67;
t26 = -pkin(2) * t46 + t68;
t25 = t31 * t59 + t32 * t63;
t24 = t31 * t63 - t32 * t59;
t16 = t24 * t58 + t25 * t62;
t15 = t24 * t62 - t25 * t58;
t12 = -pkin(5) * t24 + t17;
t4 = pkin(9) * t24 + t6;
t3 = pkin(5) * t40 - pkin(9) * t25 + t5;
t2 = t3 * t58 + t4 * t62;
t1 = t3 * t62 - t4 * t58;
t7 = m(2) * (t43 ^ 2 + t44 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t38 ^ 2) / 0.2e1 + m(4) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(6) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t22 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t12 ^ 2 + t2 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + Ifges(2,3) * t55 / 0.2e1) * t55 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t40 / 0.2e1) * t40 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t37 / 0.2e1) * t37 + (t22 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,1) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,5) * t55 + Ifges(2,1) * t48 / 0.2e1) * t48 + (-t22 * mrSges(5,1) + t14 * mrSges(5,3) + Ifges(5,4) * t32 + Ifges(5,2) * t31 / 0.2e1) * t31 + (t17 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t40 + Ifges(6,1) * t25 / 0.2e1) * t25 + (t12 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t37 + Ifges(7,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t55 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t17 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t25 + Ifges(6,6) * t40 + Ifges(6,2) * t24 / 0.2e1) * t24 + (-t12 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t37 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + t26 * mrSges(4,2) - t27 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t46) * t46 + (t38 * mrSges(3,1) + t27 * mrSges(4,1) - t28 * mrSges(4,2) - t30 * mrSges(3,3) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t41 + (Ifges(4,5) - Ifges(3,6)) * t46) * t41 + (t26 * mrSges(4,1) + t13 * mrSges(5,1) + t38 * mrSges(3,2) - t14 * mrSges(5,2) - t29 * mrSges(3,3) - t28 * mrSges(4,3) + Ifges(5,5) * t32 + Ifges(5,6) * t31 + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t42 + (-Ifges(4,4) + Ifges(3,5)) * t46 + (-Ifges(3,4) - Ifges(4,6)) * t41) * t42;
T  = t7;

% Calculate kinetic energy for
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:49
% EndTime: 2019-03-09 08:09:50
% DurationCPUTime: 0.94s
% Computational Cost: add. (2555->153), mult. (3359->209), div. (0->0), fcn. (2644->10), ass. (0->55)
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t48 = t60 * V_base(5) + t63 * V_base(4);
t54 = V_base(6) + qJD(1);
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t41 = -t48 * t59 + t54 * t62;
t42 = t48 * t62 + t54 * t59;
t56 = sin(pkin(9));
t67 = cos(pkin(9));
t32 = t56 * t41 + t42 * t67;
t47 = -t60 * V_base(4) + t63 * V_base(5);
t46 = qJD(2) - t47;
t36 = -pkin(1) * t47 - pkin(7) * t48 + V_base(3);
t52 = V_base(5) * pkin(6) + V_base(1);
t53 = -V_base(4) * pkin(6) + V_base(2);
t44 = t63 * t52 + t60 * t53;
t40 = pkin(7) * t54 + t44;
t28 = t62 * t36 - t40 * t59;
t22 = pkin(2) * t46 - qJ(3) * t42 + t28;
t29 = t59 * t36 + t62 * t40;
t25 = qJ(3) * t41 + t29;
t16 = t22 * t67 - t56 * t25;
t66 = qJD(4) - t16;
t68 = pkin(3) + qJ(5);
t10 = t32 * pkin(4) - t46 * t68 + t66;
t31 = -t41 * t67 + t42 * t56;
t43 = -t60 * t52 + t53 * t63;
t39 = -pkin(1) * t54 - t43;
t33 = -pkin(2) * t41 + qJD(3) + t39;
t65 = -qJ(4) * t32 + t33;
t13 = t31 * t68 + t65;
t55 = sin(pkin(10));
t57 = cos(pkin(10));
t6 = t55 * t10 + t57 * t13;
t17 = t56 * t22 + t67 * t25;
t15 = -t46 * qJ(4) - t17;
t5 = t57 * t10 - t13 * t55;
t11 = -pkin(4) * t31 + qJD(5) - t15;
t64 = V_base(3) ^ 2;
t61 = cos(qJ(6));
t58 = sin(qJ(6));
t30 = qJD(6) + t32;
t27 = t31 * t55 + t46 * t57;
t26 = t31 * t57 - t46 * t55;
t20 = t26 * t58 + t27 * t61;
t19 = t26 * t61 - t27 * t58;
t18 = pkin(3) * t31 + t65;
t14 = -t46 * pkin(3) + t66;
t7 = -pkin(5) * t26 + t11;
t4 = pkin(8) * t26 + t6;
t3 = pkin(5) * t32 - pkin(8) * t27 + t5;
t2 = t3 * t58 + t4 * t61;
t1 = t3 * t61 - t4 * t58;
t8 = m(2) * (t43 ^ 2 + t44 ^ 2 + t64) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t64) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t33 ^ 2) / 0.2e1 + m(3) * (t28 ^ 2 + t29 ^ 2 + t39 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t15 ^ 2 + t18 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + Ifges(2,3) * t54 / 0.2e1) * t54 + (t39 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,1) * t42 / 0.2e1) * t42 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t30 / 0.2e1) * t30 + (t11 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,5) * t54 + Ifges(2,1) * t48 / 0.2e1) * t48 + (-t39 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t42 + Ifges(3,2) * t41 / 0.2e1) * t41 + (-t11 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t30 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t54 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t30 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t33 * mrSges(4,1) + t15 * mrSges(5,1) - t18 * mrSges(5,2) - t17 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t31) * t31 + (t28 * mrSges(3,1) + t16 * mrSges(4,1) - t29 * mrSges(3,2) - t17 * mrSges(4,2) + t14 * mrSges(5,2) - t15 * mrSges(5,3) + Ifges(3,5) * t42 + Ifges(3,6) * t41 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t46 + (Ifges(5,5) - Ifges(4,6)) * t31) * t46 + (t14 * mrSges(5,1) + t5 * mrSges(6,1) + t33 * mrSges(4,2) - t6 * mrSges(6,2) - t16 * mrSges(4,3) - t18 * mrSges(5,3) + Ifges(6,5) * t27 + Ifges(6,6) * t26 + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,3) / 0.2e1) * t32 + (-Ifges(5,4) + Ifges(4,5)) * t46 + (-Ifges(4,4) - Ifges(5,6)) * t31) * t32;
T  = t8;

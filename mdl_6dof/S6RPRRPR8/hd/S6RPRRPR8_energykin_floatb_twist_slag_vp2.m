% Calculate kinetic energy for
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:30
% EndTime: 2019-03-09 05:22:31
% DurationCPUTime: 1.01s
% Computational Cost: add. (2527->153), mult. (3191->211), div. (0->0), fcn. (2396->10), ass. (0->56)
t70 = pkin(1) + pkin(7);
t62 = sin(qJ(1));
t69 = cos(qJ(1));
t47 = t62 * V_base(5) + t69 * V_base(4);
t56 = V_base(6) + qJD(1);
t52 = V_base(5) * pkin(6) + V_base(1);
t53 = -V_base(4) * pkin(6) + V_base(2);
t43 = -t62 * t52 + t69 * t53;
t67 = qJD(2) - t43;
t29 = t47 * pkin(2) - t70 * t56 + t67;
t46 = t62 * V_base(4) - t69 * V_base(5);
t68 = -qJ(2) * t47 + V_base(3);
t34 = t70 * t46 + t68;
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t22 = t61 * t29 + t65 * t34;
t45 = qJD(3) + t47;
t20 = pkin(8) * t45 + t22;
t44 = t69 * t52 + t62 * t53;
t38 = -t56 * qJ(2) - t44;
t35 = -pkin(2) * t46 - t38;
t41 = t46 * t65 - t61 * t56;
t42 = t46 * t61 + t56 * t65;
t27 = -pkin(3) * t41 - pkin(8) * t42 + t35;
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t14 = t64 * t20 + t60 * t27;
t32 = -t42 * t60 + t45 * t64;
t12 = qJ(5) * t32 + t14;
t57 = sin(pkin(10));
t58 = cos(pkin(10));
t13 = -t20 * t60 + t64 * t27;
t33 = t42 * t64 + t45 * t60;
t40 = qJD(4) - t41;
t9 = pkin(4) * t40 - qJ(5) * t33 + t13;
t6 = t58 * t12 + t57 * t9;
t5 = -t12 * t57 + t58 * t9;
t21 = t29 * t65 - t61 * t34;
t19 = -pkin(3) * t45 - t21;
t17 = -pkin(4) * t32 + qJD(5) + t19;
t66 = V_base(3) ^ 2;
t63 = cos(qJ(6));
t59 = sin(qJ(6));
t39 = qJD(6) + t40;
t37 = -t56 * pkin(1) + t67;
t36 = pkin(1) * t46 + t68;
t24 = t32 * t57 + t33 * t58;
t23 = t32 * t58 - t33 * t57;
t16 = t23 * t59 + t24 * t63;
t15 = t23 * t63 - t24 * t59;
t11 = -pkin(5) * t23 + t17;
t4 = pkin(9) * t23 + t6;
t3 = pkin(5) * t40 - pkin(9) * t24 + t5;
t2 = t3 * t59 + t4 * t63;
t1 = t3 * t63 - t4 * t59;
t7 = m(2) * (t43 ^ 2 + t44 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(4) * (t21 ^ 2 + t22 ^ 2 + t35 ^ 2) / 0.2e1 + m(3) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + m(6) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t21 * mrSges(4,1) - t22 * mrSges(4,2) + Ifges(4,3) * t45 / 0.2e1) * t45 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t39 / 0.2e1) * t39 + (t19 * mrSges(5,2) - t13 * mrSges(5,3) + Ifges(5,1) * t33 / 0.2e1) * t33 + (t17 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t35 * mrSges(4,2) - t21 * mrSges(4,3) + Ifges(4,5) * t45 + Ifges(4,1) * t42 / 0.2e1) * t42 + (-t19 * mrSges(5,1) + t14 * mrSges(5,3) + Ifges(5,4) * t33 + Ifges(5,2) * t32 / 0.2e1) * t32 + (-t17 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t24 + Ifges(6,2) * t23 / 0.2e1) * t23 + (t11 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t39 + Ifges(7,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t35 * mrSges(4,1) + t22 * mrSges(4,3) + Ifges(4,4) * t42 + Ifges(4,6) * t45 + Ifges(4,2) * t41 / 0.2e1) * t41 + (-t11 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t16 + Ifges(7,6) * t39 + Ifges(7,2) * t15 / 0.2e1) * t15 + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + t37 * mrSges(3,2) - t38 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t56) * t56 + (t37 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) - t36 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t47 + (-Ifges(3,4) + Ifges(2,5)) * t56) * t47 + (V_base(3) * mrSges(2,1) + t38 * mrSges(3,1) - t36 * mrSges(3,2) - t44 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t46 + (Ifges(3,5) - Ifges(2,6)) * t56 + (-Ifges(2,4) - Ifges(3,6)) * t47) * t46 + (t13 * mrSges(5,1) + t5 * mrSges(6,1) - t14 * mrSges(5,2) - t6 * mrSges(6,2) + Ifges(5,5) * t33 + Ifges(6,5) * t24 + Ifges(5,6) * t32 + Ifges(6,6) * t23 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t40) * t40;
T  = t7;

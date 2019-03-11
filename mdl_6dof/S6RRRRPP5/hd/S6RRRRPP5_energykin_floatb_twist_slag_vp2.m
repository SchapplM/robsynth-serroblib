% Calculate kinetic energy for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:00
% EndTime: 2019-03-09 21:05:01
% DurationCPUTime: 1.01s
% Computational Cost: add. (2311->149), mult. (2927->196), div. (0->0), fcn. (2268->8), ass. (0->50)
t67 = -pkin(4) - pkin(5);
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t45 = -t57 * V_base(4) + t60 * V_base(5);
t46 = t57 * V_base(5) + t60 * V_base(4);
t32 = -pkin(1) * t45 - pkin(7) * t46 + V_base(3);
t51 = pkin(6) * V_base(5) + V_base(1);
t52 = -pkin(6) * V_base(4) + V_base(2);
t43 = t51 * t60 + t52 * t57;
t53 = V_base(6) + qJD(1);
t38 = pkin(7) * t53 + t43;
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t26 = t32 * t56 + t38 * t59;
t44 = qJD(2) - t45;
t23 = pkin(8) * t44 + t26;
t42 = -t51 * t57 + t52 * t60;
t37 = -pkin(1) * t53 - t42;
t40 = -t46 * t56 + t53 * t59;
t41 = t46 * t59 + t53 * t56;
t24 = -pkin(2) * t40 - pkin(8) * t41 + t37;
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t14 = -t23 * t55 + t24 * t58;
t29 = t41 * t58 + t44 * t55;
t39 = qJD(3) - t40;
t10 = pkin(3) * t39 - pkin(9) * t29 + t14;
t15 = t23 * t58 + t24 * t55;
t28 = -t41 * t55 + t44 * t58;
t13 = pkin(9) * t28 + t15;
t54 = sin(qJ(4));
t66 = cos(qJ(4));
t6 = t10 * t54 + t13 * t66;
t25 = t32 * t59 - t38 * t56;
t36 = qJD(4) + t39;
t4 = qJ(5) * t36 + t6;
t65 = pkin(2) * t44 + t25;
t5 = t10 * t66 - t13 * t54;
t64 = qJD(5) - t5;
t63 = pkin(3) * t28 + t65;
t18 = t28 * t54 + t29 * t66;
t62 = qJ(5) * t18 + t63;
t61 = V_base(3) ^ 2;
t17 = -t28 * t66 + t29 * t54;
t8 = pkin(4) * t17 - t62;
t7 = t17 * t67 + qJD(6) + t62;
t3 = -pkin(4) * t36 + t64;
t2 = qJ(6) * t17 + t4;
t1 = -qJ(6) * t18 + t67 * t36 + t64;
t9 = m(2) * (t42 ^ 2 + t43 ^ 2 + t61) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t61) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t37 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t65 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t63 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t42 * mrSges(2,1) - t43 * mrSges(2,2) + Ifges(2,3) * t53 / 0.2e1) * t53 + (t25 * mrSges(3,1) - t26 * mrSges(3,2) + Ifges(3,3) * t44 / 0.2e1) * t44 + (t14 * mrSges(4,1) - t15 * mrSges(4,2) + Ifges(4,3) * t39 / 0.2e1) * t39 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t42 * mrSges(2,3) + Ifges(2,5) * t53 + Ifges(2,1) * t46 / 0.2e1) * t46 + (t37 * mrSges(3,2) - t25 * mrSges(3,3) + t44 * Ifges(3,5) + Ifges(3,1) * t41 / 0.2e1) * t41 + (-t65 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,5) * t39 + Ifges(4,1) * t29 / 0.2e1) * t29 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t43 * mrSges(2,3) + Ifges(2,4) * t46 + Ifges(2,6) * t53 + Ifges(2,2) * t45 / 0.2e1) * t45 + (-t37 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t41 + Ifges(3,6) * t44 + Ifges(3,2) * t40 / 0.2e1) * t40 + (t65 * mrSges(4,1) + t15 * mrSges(4,3) + Ifges(4,4) * t29 + Ifges(4,6) * t39 + Ifges(4,2) * t28 / 0.2e1) * t28 + (t5 * mrSges(5,1) - t3 * mrSges(6,1) - t1 * mrSges(7,1) - t6 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t36) * t36 + (-t63 * mrSges(5,2) + t3 * mrSges(6,2) + t7 * mrSges(7,2) - t5 * mrSges(5,3) - t8 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t18 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t36) * t18 + (-t63 * mrSges(5,1) + t8 * mrSges(6,1) - t7 * mrSges(7,1) - t4 * mrSges(6,2) - t6 * mrSges(5,3) + t2 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t17 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t36 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t18) * t17;
T  = t9;

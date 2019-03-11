% Calculate kinetic energy for
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:50
% EndTime: 2019-03-09 15:24:51
% DurationCPUTime: 1.03s
% Computational Cost: add. (2973->153), mult. (3927->211), div. (0->0), fcn. (3148->10), ass. (0->56)
t69 = pkin(4) + pkin(9);
t60 = sin(qJ(1));
t64 = cos(qJ(1));
t47 = -t60 * V_base(4) + t64 * V_base(5);
t48 = t60 * V_base(5) + t64 * V_base(4);
t37 = -pkin(1) * t47 - pkin(7) * t48 + V_base(3);
t52 = pkin(6) * V_base(5) + V_base(1);
t53 = -pkin(6) * V_base(4) + V_base(2);
t44 = t52 * t64 + t53 * t60;
t55 = V_base(6) + qJD(1);
t40 = pkin(7) * t55 + t44;
t59 = sin(qJ(2));
t63 = cos(qJ(2));
t29 = t37 * t63 - t40 * t59;
t42 = t48 * t63 + t55 * t59;
t46 = qJD(2) - t47;
t26 = pkin(2) * t46 - pkin(8) * t42 + t29;
t30 = t37 * t59 + t40 * t63;
t41 = -t48 * t59 + t55 * t63;
t28 = pkin(8) * t41 + t30;
t58 = sin(qJ(3));
t62 = cos(qJ(3));
t16 = t26 * t62 - t28 * t58;
t33 = t41 * t58 + t42 * t62;
t45 = qJD(3) + t46;
t12 = pkin(3) * t45 - qJ(4) * t33 + t16;
t17 = t26 * t58 + t28 * t62;
t32 = t41 * t62 - t42 * t58;
t15 = qJ(4) * t32 + t17;
t56 = sin(pkin(10));
t68 = cos(pkin(10));
t8 = t12 * t56 + t15 * t68;
t43 = -t52 * t60 + t53 * t64;
t6 = -qJ(5) * t45 - t8;
t7 = t12 * t68 - t15 * t56;
t39 = -pkin(1) * t55 - t43;
t67 = qJD(5) - t7;
t22 = t32 * t56 + t33 * t68;
t34 = -pkin(2) * t41 + t39;
t23 = -pkin(3) * t32 + qJD(4) + t34;
t66 = -qJ(5) * t22 + t23;
t65 = V_base(3) ^ 2;
t61 = cos(qJ(6));
t57 = sin(qJ(6));
t21 = -t32 * t68 + t33 * t56;
t20 = qJD(6) + t22;
t19 = t21 * t57 + t45 * t61;
t18 = t21 * t61 - t45 * t57;
t10 = pkin(4) * t21 + t66;
t9 = t21 * t69 + t66;
t5 = -pkin(4) * t45 + t67;
t4 = -pkin(5) * t21 - t6;
t3 = pkin(5) * t22 - t45 * t69 + t67;
t2 = t3 * t57 + t61 * t9;
t1 = t3 * t61 - t57 * t9;
t11 = m(2) * (t43 ^ 2 + t44 ^ 2 + t65) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t65) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t39 ^ 2) / 0.2e1 + m(5) * (t23 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t34 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + Ifges(2,3) * t55 / 0.2e1) * t55 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t34 * mrSges(4,2) - t16 * mrSges(4,3) + Ifges(4,1) * t33 / 0.2e1) * t33 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,5) * t55 + Ifges(2,1) * t48 / 0.2e1) * t48 + (t39 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t42 / 0.2e1) * t42 + (-t34 * mrSges(4,1) + t17 * mrSges(4,3) + Ifges(4,4) * t33 + Ifges(4,2) * t32 / 0.2e1) * t32 + (t4 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t20 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t55 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t39 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t42 + Ifges(3,6) * t46 + Ifges(3,2) * t41 / 0.2e1) * t41 + (-t4 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t20 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t5 * mrSges(6,1) + t23 * mrSges(5,2) - t7 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t22) * t22 + (t23 * mrSges(5,1) + t6 * mrSges(6,1) - t10 * mrSges(6,2) - t8 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t21 + (-Ifges(5,4) - Ifges(6,6)) * t22) * t21 + (t16 * mrSges(4,1) + t7 * mrSges(5,1) - t17 * mrSges(4,2) - t8 * mrSges(5,2) + t5 * mrSges(6,2) - t6 * mrSges(6,3) + Ifges(4,5) * t33 + Ifges(4,6) * t32 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t45 + (-Ifges(6,4) + Ifges(5,5)) * t22 + (Ifges(6,5) - Ifges(5,6)) * t21) * t45;
T  = t11;

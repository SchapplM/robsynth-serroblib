% Calculate kinetic energy for
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:19
% EndTime: 2019-03-09 03:44:20
% DurationCPUTime: 1.03s
% Computational Cost: add. (2425->153), mult. (3491->211), div. (0->0), fcn. (2756->10), ass. (0->56)
t69 = pkin(3) + pkin(8);
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t47 = -t61 * V_base(4) + t65 * V_base(5);
t48 = t61 * V_base(5) + t65 * V_base(4);
t56 = sin(pkin(10));
t57 = cos(pkin(10));
t43 = t47 * t56 + t48 * t57;
t55 = V_base(6) + qJD(1);
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t35 = t43 * t64 + t55 * t60;
t42 = t47 * t57 - t48 * t56;
t41 = qJD(3) - t42;
t53 = pkin(6) * V_base(5) + V_base(1);
t54 = -pkin(6) * V_base(4) + V_base(2);
t44 = -t53 * t61 + t54 * t65;
t36 = pkin(1) * t55 - qJ(2) * t48 + t44;
t45 = t53 * t65 + t54 * t61;
t40 = qJ(2) * t47 + t45;
t30 = t36 * t56 + t40 * t57;
t24 = pkin(7) * t55 + t30;
t46 = -pkin(1) * t47 + qJD(2) + V_base(3);
t26 = -pkin(2) * t42 - pkin(7) * t43 + t46;
t16 = -t24 * t60 + t26 * t64;
t68 = qJD(4) - t16;
t10 = pkin(4) * t35 - t41 * t69 + t68;
t34 = t43 * t60 - t55 * t64;
t29 = t36 * t57 - t40 * t56;
t23 = -pkin(2) * t55 - t29;
t67 = -qJ(4) * t35 + t23;
t13 = t34 * t69 + t67;
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t6 = t10 * t59 + t13 * t63;
t17 = t24 * t64 + t26 * t60;
t15 = -qJ(4) * t41 - t17;
t5 = t10 * t63 - t13 * t59;
t11 = -pkin(4) * t34 - t15;
t32 = qJD(5) + t35;
t66 = V_base(3) ^ 2;
t62 = cos(qJ(6));
t58 = sin(qJ(6));
t31 = qJD(6) + t32;
t28 = t34 * t59 + t41 * t63;
t27 = t34 * t63 - t41 * t59;
t20 = t27 * t58 + t28 * t62;
t19 = t27 * t62 - t28 * t58;
t18 = pkin(3) * t34 + t67;
t14 = -pkin(3) * t41 + t68;
t7 = -pkin(5) * t27 + t11;
t4 = pkin(9) * t27 + t6;
t3 = pkin(5) * t32 - pkin(9) * t28 + t5;
t2 = t3 * t58 + t4 * t62;
t1 = t3 * t62 - t4 * t58;
t8 = m(2) * (t44 ^ 2 + t45 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t23 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t15 ^ 2 + t18 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,1) * t48 / 0.2e1) * t48 + (t46 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,1) * t43 / 0.2e1) * t43 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t32 / 0.2e1) * t32 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t46 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t43 + Ifges(3,2) * t42 / 0.2e1) * t42 + (t11 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t32 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t31 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t11 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t32 + Ifges(6,2) * t27 / 0.2e1) * t27 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t31 + Ifges(7,2) * t19 / 0.2e1) * t19 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + t14 * mrSges(5,2) - t15 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t41) * t41 + (t14 * mrSges(5,1) + t23 * mrSges(4,2) - t16 * mrSges(4,3) - t18 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t35 + (-Ifges(5,4) + Ifges(4,5)) * t41) * t35 + (t44 * mrSges(2,1) + t29 * mrSges(3,1) - t45 * mrSges(2,2) - t30 * mrSges(3,2) + Ifges(2,5) * t48 + Ifges(3,5) * t43 + Ifges(2,6) * t47 + Ifges(3,6) * t42 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * t55) * t55 + (t23 * mrSges(4,1) + t15 * mrSges(5,1) - t18 * mrSges(5,2) - t17 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t34 + (Ifges(5,5) - Ifges(4,6)) * t41 + (-Ifges(4,4) - Ifges(5,6)) * t35) * t34;
T  = t8;

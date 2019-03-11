% Calculate kinetic energy for
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:02
% EndTime: 2019-03-09 02:50:03
% DurationCPUTime: 0.95s
% Computational Cost: add. (2527->153), mult. (3359->209), div. (0->0), fcn. (2644->10), ass. (0->55)
t62 = sin(qJ(1));
t64 = cos(qJ(1));
t48 = t62 * V_base(5) + t64 * V_base(4);
t55 = V_base(6) + qJD(1);
t57 = sin(pkin(9));
t59 = cos(pkin(9));
t41 = -t48 * t57 + t55 * t59;
t42 = t48 * t59 + t55 * t57;
t61 = sin(qJ(3));
t69 = cos(qJ(3));
t32 = t61 * t41 + t42 * t69;
t47 = t62 * V_base(4) - t64 * V_base(5);
t46 = qJD(3) + t47;
t36 = pkin(1) * t47 - qJ(2) * t48 + V_base(3);
t52 = V_base(5) * pkin(6) + V_base(1);
t53 = -V_base(4) * pkin(6) + V_base(2);
t44 = t64 * t52 + t62 * t53;
t40 = qJ(2) * t55 + t44;
t28 = t59 * t36 - t40 * t57;
t22 = pkin(2) * t47 - pkin(7) * t42 + t28;
t29 = t57 * t36 + t59 * t40;
t25 = pkin(7) * t41 + t29;
t16 = t22 * t69 - t61 * t25;
t67 = qJD(4) - t16;
t68 = pkin(3) + qJ(5);
t10 = t32 * pkin(4) - t46 * t68 + t67;
t31 = -t41 * t69 + t42 * t61;
t43 = -t62 * t52 + t53 * t64;
t38 = -pkin(1) * t55 + qJD(2) - t43;
t33 = -pkin(2) * t41 + t38;
t66 = -qJ(4) * t32 + t33;
t13 = t31 * t68 + t66;
t56 = sin(pkin(10));
t58 = cos(pkin(10));
t6 = t56 * t10 + t58 * t13;
t17 = t61 * t22 + t69 * t25;
t15 = -t46 * qJ(4) - t17;
t5 = t58 * t10 - t13 * t56;
t11 = -pkin(4) * t31 + qJD(5) - t15;
t65 = V_base(3) ^ 2;
t63 = cos(qJ(6));
t60 = sin(qJ(6));
t30 = qJD(6) + t32;
t27 = t31 * t56 + t46 * t58;
t26 = t31 * t58 - t46 * t56;
t20 = t26 * t60 + t27 * t63;
t19 = t26 * t63 - t27 * t60;
t18 = pkin(3) * t31 + t66;
t14 = -t46 * pkin(3) + t67;
t7 = -pkin(5) * t26 + t11;
t4 = pkin(8) * t26 + t6;
t3 = pkin(5) * t32 - pkin(8) * t27 + t5;
t2 = t3 * t60 + t4 * t63;
t1 = t3 * t63 - t4 * t60;
t8 = m(2) * (t43 ^ 2 + t44 ^ 2 + t65) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t65) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t33 ^ 2) / 0.2e1 + m(3) * (t28 ^ 2 + t29 ^ 2 + t38 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t15 ^ 2 + t18 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + Ifges(2,3) * t55 / 0.2e1) * t55 + (t38 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,1) * t42 / 0.2e1) * t42 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t30 / 0.2e1) * t30 + (t11 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,5) * t55 + Ifges(2,1) * t48 / 0.2e1) * t48 + (-t38 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t42 + Ifges(3,2) * t41 / 0.2e1) * t41 + (-t11 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t27 + Ifges(6,2) * t26 / 0.2e1) * t26 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t30 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + t20 * Ifges(7,4) + t30 * Ifges(7,6) + Ifges(7,2) * t19 / 0.2e1) * t19 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + t14 * mrSges(5,2) - t15 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t46) * t46 + (t33 * mrSges(4,1) + t15 * mrSges(5,1) - t18 * mrSges(5,2) - t17 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t31 + (Ifges(5,5) - Ifges(4,6)) * t46) * t31 + (V_base(3) * mrSges(2,1) + t28 * mrSges(3,1) - t29 * mrSges(3,2) - t44 * mrSges(2,3) - Ifges(2,4) * t48 + Ifges(3,5) * t42 - Ifges(2,6) * t55 + Ifges(3,6) * t41 + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t47) * t47 + (t14 * mrSges(5,1) + t5 * mrSges(6,1) + t33 * mrSges(4,2) - t6 * mrSges(6,2) - t16 * mrSges(4,3) - t18 * mrSges(5,3) + Ifges(6,5) * t27 + Ifges(6,6) * t26 + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,3) / 0.2e1) * t32 + (-Ifges(5,4) + Ifges(4,5)) * t46 + (-Ifges(4,4) - Ifges(5,6)) * t31) * t32;
T  = t8;

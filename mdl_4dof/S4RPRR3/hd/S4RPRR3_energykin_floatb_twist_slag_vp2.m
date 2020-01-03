% Calculate kinetic energy for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:03
% EndTime: 2019-12-31 16:49:03
% DurationCPUTime: 0.56s
% Computational Cost: add. (933->108), mult. (1357->158), div. (0->0), fcn. (996->8), ass. (0->40)
t38 = V_base(5) * pkin(4) + V_base(1);
t39 = -V_base(4) * pkin(4) + V_base(2);
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t29 = -t38 * t45 + t48 * t39;
t34 = t45 * V_base(5) + t48 * V_base(4);
t40 = V_base(6) + qJD(1);
t21 = pkin(1) * t40 - qJ(2) * t34 + t29;
t30 = t48 * t38 + t45 * t39;
t33 = -t45 * V_base(4) + t48 * V_base(5);
t25 = qJ(2) * t33 + t30;
t41 = sin(pkin(7));
t42 = cos(pkin(7));
t17 = t41 * t21 + t42 * t25;
t12 = pkin(5) * t40 + t17;
t27 = t33 * t42 - t41 * t34;
t28 = t33 * t41 + t34 * t42;
t31 = -pkin(1) * t33 + qJD(2) + V_base(3);
t15 = -pkin(2) * t27 - pkin(5) * t28 + t31;
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t6 = t47 * t12 + t44 * t15;
t5 = -t12 * t44 + t47 * t15;
t16 = t21 * t42 - t41 * t25;
t26 = qJD(3) - t27;
t11 = -pkin(2) * t40 - t16;
t49 = V_base(3) ^ 2;
t46 = cos(qJ(4));
t43 = sin(qJ(4));
t24 = qJD(4) + t26;
t20 = t28 * t47 + t40 * t44;
t19 = -t28 * t44 + t40 * t47;
t9 = t19 * t43 + t20 * t46;
t8 = t19 * t46 - t20 * t43;
t7 = -pkin(3) * t19 + t11;
t4 = pkin(6) * t19 + t6;
t3 = pkin(3) * t26 - pkin(6) * t20 + t5;
t2 = t3 * t43 + t4 * t46;
t1 = t3 * t46 - t4 * t43;
t10 = m(2) * (t29 ^ 2 + t30 ^ 2 + t49) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t49) / 0.2e1 + m(3) * (t16 ^ 2 + t17 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t7 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,1) * t9 / 0.2e1) * t9 + (V_base(3) * mrSges(2,2) - t29 * mrSges(2,3) + Ifges(2,1) * t34 / 0.2e1) * t34 + (t31 * mrSges(3,2) - t16 * mrSges(3,3) + Ifges(3,1) * t28 / 0.2e1) * t28 + (t5 * mrSges(4,1) - t6 * mrSges(4,2) + Ifges(4,3) * t26 / 0.2e1) * t26 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t7 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,2) * t8 / 0.2e1) * t8 + (-V_base(3) * mrSges(2,1) + t30 * mrSges(2,3) + Ifges(2,4) * t34 + Ifges(2,2) * t33 / 0.2e1) * t33 + (-t31 * mrSges(3,1) + t17 * mrSges(3,3) + Ifges(3,4) * t28 + Ifges(3,2) * t27 / 0.2e1) * t27 + (t11 * mrSges(4,2) - t5 * mrSges(4,3) + Ifges(4,5) * t26 + Ifges(4,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t8 + Ifges(5,3) * t24 / 0.2e1) * t24 + (-t11 * mrSges(4,1) + t6 * mrSges(4,3) + Ifges(4,4) * t20 + Ifges(4,6) * t26 + Ifges(4,2) * t19 / 0.2e1) * t19 + (t29 * mrSges(2,1) + t16 * mrSges(3,1) - t30 * mrSges(2,2) - t17 * mrSges(3,2) + Ifges(2,5) * t34 + Ifges(3,5) * t28 + Ifges(2,6) * t33 + Ifges(3,6) * t27 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t40) * t40;
T = t10;

% Calculate kinetic energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:43
% EndTime: 2019-12-31 17:40:43
% DurationCPUTime: 0.70s
% Computational Cost: add. (1017->125), mult. (1504->160), div. (0->0), fcn. (1088->6), ass. (0->39)
t52 = -pkin(3) - pkin(4);
t51 = cos(qJ(2));
t50 = cos(qJ(3));
t39 = V_base(5) * qJ(1) + V_base(1);
t40 = -V_base(4) * qJ(1) + V_base(2);
t43 = sin(pkin(7));
t44 = cos(pkin(7));
t29 = -t39 * t43 + t44 * t40;
t34 = t43 * V_base(5) + t44 * V_base(4);
t21 = V_base(6) * pkin(1) - pkin(5) * t34 + t29;
t30 = t44 * t39 + t43 * t40;
t33 = -t43 * V_base(4) + t44 * V_base(5);
t24 = pkin(5) * t33 + t30;
t46 = sin(qJ(2));
t16 = t46 * t21 + t51 * t24;
t41 = V_base(6) + qJD(2);
t13 = pkin(6) * t41 + t16;
t27 = t51 * t33 - t34 * t46;
t28 = t46 * t33 + t51 * t34;
t42 = V_base(3) + qJD(1);
t31 = -pkin(1) * t33 + t42;
t14 = -pkin(2) * t27 - pkin(6) * t28 + t31;
t45 = sin(qJ(3));
t7 = t50 * t13 + t45 * t14;
t15 = t51 * t21 - t46 * t24;
t26 = qJD(3) - t27;
t5 = t26 * qJ(4) + t7;
t49 = pkin(2) * t41 + t15;
t6 = -t45 * t13 + t50 * t14;
t48 = qJD(4) - t6;
t18 = t50 * t28 + t45 * t41;
t47 = qJ(4) * t18 + t49;
t17 = t28 * t45 - t50 * t41;
t8 = pkin(3) * t17 - t47;
t4 = -t26 * pkin(3) + t48;
t3 = t52 * t17 + qJD(5) + t47;
t2 = qJ(5) * t17 + t5;
t1 = -t18 * qJ(5) + t52 * t26 + t48;
t9 = m(2) * (t29 ^ 2 + t30 ^ 2 + t42 ^ 2) / 0.2e1 + m(3) * (t15 ^ 2 + t16 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t49 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t15 * mrSges(3,1) - t16 * mrSges(3,2) + Ifges(3,3) * t41 / 0.2e1) * t41 + (t42 * mrSges(2,2) - t29 * mrSges(2,3) + Ifges(2,1) * t34 / 0.2e1) * t34 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t42 * mrSges(2,1) + t30 * mrSges(2,3) + Ifges(2,4) * t34 + Ifges(2,2) * t33 / 0.2e1) * t33 + (t31 * mrSges(3,2) - t15 * mrSges(3,3) + Ifges(3,5) * t41 + Ifges(3,1) * t28 / 0.2e1) * t28 + (-t31 * mrSges(3,1) + t16 * mrSges(3,3) + Ifges(3,4) * t28 + Ifges(3,6) * t41 + Ifges(3,2) * t27 / 0.2e1) * t27 + (V_base(2) * mrSges(1,1) + t29 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t30 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t34 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t33 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t6 * mrSges(4,1) - t4 * mrSges(5,1) - t1 * mrSges(6,1) - t7 * mrSges(4,2) + t2 * mrSges(6,2) + t5 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t26) * t26 + (-t49 * mrSges(4,2) + t4 * mrSges(5,2) + t3 * mrSges(6,2) - t6 * mrSges(4,3) - t8 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t18 + (Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t26) * t18 + (-t49 * mrSges(4,1) + t8 * mrSges(5,1) - t3 * mrSges(6,1) - t5 * mrSges(5,2) - t7 * mrSges(4,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t17 + (-Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t26 + (-Ifges(4,4) + Ifges(6,4) + Ifges(5,5)) * t18) * t17;
T = t9;

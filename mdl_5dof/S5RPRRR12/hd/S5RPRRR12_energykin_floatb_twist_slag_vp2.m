% Calculate kinetic energy for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:25
% EndTime: 2019-12-31 19:12:25
% DurationCPUTime: 0.78s
% Computational Cost: add. (1215->129), mult. (1528->178), div. (0->0), fcn. (1056->8), ass. (0->47)
t56 = pkin(1) + pkin(6);
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t36 = t48 * V_base(5) + t52 * V_base(4);
t44 = V_base(6) + qJD(1);
t40 = V_base(5) * pkin(5) + V_base(1);
t41 = -V_base(4) * pkin(5) + V_base(2);
t31 = -t48 * t40 + t41 * t52;
t54 = qJD(2) - t31;
t22 = pkin(2) * t36 - t44 * t56 + t54;
t35 = t48 * V_base(4) - t52 * V_base(5);
t55 = -qJ(2) * t36 + V_base(3);
t24 = t35 * t56 + t55;
t47 = sin(qJ(3));
t51 = cos(qJ(3));
t13 = t47 * t22 + t51 * t24;
t29 = t35 * t51 - t44 * t47;
t11 = pkin(7) * t29 + t13;
t46 = sin(qJ(4));
t50 = cos(qJ(4));
t12 = t51 * t22 - t24 * t47;
t30 = t35 * t47 + t44 * t51;
t34 = qJD(3) + t36;
t9 = pkin(3) * t34 - pkin(7) * t30 + t12;
t6 = t50 * t11 + t46 * t9;
t32 = t52 * t40 + t48 * t41;
t28 = -t44 * qJ(2) - t32;
t5 = -t11 * t46 + t50 * t9;
t18 = t29 * t50 - t30 * t46;
t25 = -pkin(2) * t35 - t28;
t16 = -pkin(3) * t29 + t25;
t53 = V_base(3) ^ 2;
t49 = cos(qJ(5));
t45 = sin(qJ(5));
t33 = qJD(4) + t34;
t27 = -pkin(1) * t44 + t54;
t26 = pkin(1) * t35 + t55;
t19 = t29 * t46 + t30 * t50;
t17 = qJD(5) - t18;
t15 = t19 * t49 + t33 * t45;
t14 = -t19 * t45 + t33 * t49;
t7 = -pkin(4) * t18 - pkin(8) * t19 + t16;
t4 = pkin(8) * t33 + t6;
t3 = -pkin(4) * t33 - t5;
t2 = t4 * t49 + t45 * t7;
t1 = -t4 * t45 + t49 * t7;
t8 = m(2) * (t31 ^ 2 + t32 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t25 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(5) * (t16 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t12 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,3) * t34 / 0.2e1) * t34 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t33 / 0.2e1) * t33 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t17 / 0.2e1) * t17 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t25 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,5) * t34 + Ifges(4,1) * t30 / 0.2e1) * t30 + (t16 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t33 + Ifges(5,1) * t19 / 0.2e1) * t19 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t17 + Ifges(6,1) * t15 / 0.2e1) * t15 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t25 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,6) * t34 + Ifges(4,2) * t29 / 0.2e1) * t29 + (-t16 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t19 + Ifges(5,6) * t33 + Ifges(5,2) * t18 / 0.2e1) * t18 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t17 + Ifges(6,2) * t14 / 0.2e1) * t14 + (t31 * mrSges(2,1) - t32 * mrSges(2,2) + t27 * mrSges(3,2) - t28 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t44) * t44 + (t27 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t31 * mrSges(2,3) - t26 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t36 + (-Ifges(3,4) + Ifges(2,5)) * t44) * t36 + (V_base(3) * mrSges(2,1) + t28 * mrSges(3,1) - t26 * mrSges(3,2) - t32 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t35 + (Ifges(3,5) - Ifges(2,6)) * t44 + (-Ifges(2,4) - Ifges(3,6)) * t36) * t35;
T = t8;

% Calculate kinetic energy for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:28:56
% EndTime: 2019-12-31 20:28:57
% DurationCPUTime: 0.78s
% Computational Cost: add. (1229->129), mult. (1554->178), div. (0->0), fcn. (1112->8), ass. (0->45)
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t38 = -t53 * V_base(4) + t56 * V_base(5);
t39 = t53 * V_base(5) + t56 * V_base(4);
t24 = -pkin(1) * t38 - pkin(6) * t39 + V_base(3);
t45 = V_base(5) * pkin(5) + V_base(1);
t46 = -V_base(4) * pkin(5) + V_base(2);
t34 = t56 * t45 + t53 * t46;
t49 = V_base(6) + qJD(1);
t28 = pkin(6) * t49 + t34;
t52 = sin(qJ(2));
t59 = cos(qJ(2));
t19 = t52 * t24 + t59 * t28;
t37 = qJD(2) - t38;
t14 = t37 * qJ(3) + t19;
t31 = t39 * t52 - t49 * t59;
t11 = pkin(7) * t31 + t14;
t51 = sin(qJ(4));
t55 = cos(qJ(4));
t32 = t39 * t59 + t52 * t49;
t18 = t24 * t59 - t52 * t28;
t58 = qJD(3) - t18;
t9 = -t32 * pkin(7) + (-pkin(2) - pkin(3)) * t37 + t58;
t6 = t55 * t11 + t51 * t9;
t33 = -t53 * t45 + t56 * t46;
t27 = -t49 * pkin(1) - t33;
t15 = t31 * pkin(2) - t32 * qJ(3) + t27;
t5 = -t11 * t51 + t55 * t9;
t21 = t31 * t55 - t32 * t51;
t12 = -pkin(3) * t31 - t15;
t57 = V_base(3) ^ 2;
t54 = cos(qJ(5));
t50 = sin(qJ(5));
t36 = qJD(4) - t37;
t22 = t31 * t51 + t32 * t55;
t20 = qJD(5) - t21;
t17 = t22 * t54 + t36 * t50;
t16 = -t22 * t50 + t36 * t54;
t13 = -t37 * pkin(2) + t58;
t7 = -pkin(4) * t21 - pkin(8) * t22 + t12;
t4 = pkin(8) * t36 + t6;
t3 = -pkin(4) * t36 - t5;
t2 = t4 * t54 + t50 * t7;
t1 = -t4 * t50 + t54 * t7;
t8 = m(2) * (t33 ^ 2 + t34 ^ 2 + t57) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t57) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t15 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t49 / 0.2e1) * t49 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t36 / 0.2e1) * t36 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,5) * t49 + Ifges(2,1) * t39 / 0.2e1) * t39 + (t12 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t36 + Ifges(5,1) * t22 / 0.2e1) * t22 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t20 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t39 + Ifges(2,6) * t49 + Ifges(2,2) * t38 / 0.2e1) * t38 + (-t12 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t22 + Ifges(5,6) * t36 + Ifges(5,2) * t21 / 0.2e1) * t21 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t20 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t18 * mrSges(3,1) - t13 * mrSges(4,1) - t19 * mrSges(3,2) + t14 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t37) * t37 + (t27 * mrSges(3,2) + t13 * mrSges(4,2) - t18 * mrSges(3,3) - t15 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t32 + (Ifges(4,4) + Ifges(3,5)) * t37) * t32 + (t27 * mrSges(3,1) + t15 * mrSges(4,1) - t14 * mrSges(4,2) - t19 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t31 + (-Ifges(3,6) + Ifges(4,6)) * t37 + (-Ifges(3,4) + Ifges(4,5)) * t32) * t31;
T = t8;

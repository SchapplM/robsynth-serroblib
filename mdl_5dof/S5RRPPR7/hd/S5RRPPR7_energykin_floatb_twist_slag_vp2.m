% Calculate kinetic energy for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:51
% EndTime: 2019-12-31 19:34:51
% DurationCPUTime: 0.71s
% Computational Cost: add. (1351->129), mult. (1774->176), div. (0->0), fcn. (1316->8), ass. (0->46)
t55 = pkin(3) + pkin(7);
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t36 = -t47 * V_base(4) + t50 * V_base(5);
t37 = t47 * V_base(5) + t50 * V_base(4);
t26 = -pkin(1) * t36 - pkin(6) * t37 + V_base(3);
t41 = pkin(5) * V_base(5) + V_base(1);
t42 = -pkin(5) * V_base(4) + V_base(2);
t34 = t50 * t41 + t47 * t42;
t43 = V_base(6) + qJD(1);
t30 = pkin(6) * t43 + t34;
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t18 = t49 * t26 - t30 * t46;
t32 = t37 * t49 + t43 * t46;
t35 = qJD(2) - t36;
t12 = pkin(2) * t35 - qJ(3) * t32 + t18;
t19 = t46 * t26 + t49 * t30;
t31 = -t37 * t46 + t43 * t49;
t15 = qJ(3) * t31 + t19;
t44 = sin(pkin(8));
t54 = cos(pkin(8));
t9 = t44 * t12 + t54 * t15;
t33 = -t47 * t41 + t42 * t50;
t7 = -qJ(4) * t35 - t9;
t8 = t12 * t54 - t44 * t15;
t29 = -pkin(1) * t43 - t33;
t53 = qJD(4) - t8;
t22 = t44 * t31 + t32 * t54;
t23 = -pkin(2) * t31 + qJD(3) + t29;
t52 = -qJ(4) * t22 + t23;
t51 = V_base(3) ^ 2;
t48 = cos(qJ(5));
t45 = sin(qJ(5));
t21 = -t31 * t54 + t32 * t44;
t20 = qJD(5) + t22;
t17 = t21 * t45 + t35 * t48;
t16 = t21 * t48 - t35 * t45;
t10 = pkin(3) * t21 + t52;
t6 = -t35 * pkin(3) + t53;
t5 = t21 * t55 + t52;
t4 = -pkin(4) * t21 - t7;
t3 = t22 * pkin(4) - t35 * t55 + t53;
t2 = t3 * t45 + t48 * t5;
t1 = t3 * t48 - t45 * t5;
t11 = m(2) * (t33 ^ 2 + t34 ^ 2 + t51) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t51) / 0.2e1 + m(4) * (t23 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t29 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t43 / 0.2e1) * t43 + (t29 * mrSges(3,2) - t18 * mrSges(3,3) + Ifges(3,1) * t32 / 0.2e1) * t32 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,5) * t43 + Ifges(2,1) * t37 / 0.2e1) * t37 + (-t29 * mrSges(3,1) + t19 * mrSges(3,3) + Ifges(3,4) * t32 + Ifges(3,2) * t31 / 0.2e1) * t31 + (t4 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t20 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t37 + Ifges(2,6) * t43 + Ifges(2,2) * t36 / 0.2e1) * t36 + (-t4 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t20 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t6 * mrSges(5,1) + t23 * mrSges(4,2) - t8 * mrSges(4,3) - t10 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t22) * t22 + (t23 * mrSges(4,1) + t7 * mrSges(5,1) - t10 * mrSges(5,2) - t9 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t21 + (-Ifges(4,4) - Ifges(5,6)) * t22) * t21 + (t18 * mrSges(3,1) + t8 * mrSges(4,1) - t19 * mrSges(3,2) - t9 * mrSges(4,2) + t6 * mrSges(5,2) - t7 * mrSges(5,3) + Ifges(3,5) * t32 + Ifges(3,6) * t31 + (Ifges(3,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t35 + (-Ifges(5,4) + Ifges(4,5)) * t22 + (Ifges(5,5) - Ifges(4,6)) * t21) * t35;
T = t11;

% Calculate kinetic energy for
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:37
% EndTime: 2019-12-05 16:08:38
% DurationCPUTime: 0.68s
% Computational Cost: add. (1383->128), mult. (2024->175), div. (0->0), fcn. (1524->8), ass. (0->42)
t46 = sin(pkin(8));
t53 = cos(pkin(8));
t47 = sin(pkin(7));
t48 = cos(pkin(7));
t38 = -t47 * V_base(4) + t48 * V_base(5);
t39 = t47 * V_base(5) + t48 * V_base(4);
t45 = V_base(3) + qJD(1);
t27 = -pkin(1) * t38 - pkin(5) * t39 + t45;
t43 = V_base(5) * qJ(1) + V_base(1);
t44 = -V_base(4) * qJ(1) + V_base(2);
t36 = t48 * t43 + t47 * t44;
t31 = V_base(6) * pkin(5) + t36;
t50 = sin(qJ(2));
t52 = cos(qJ(2));
t22 = t50 * t27 + t52 * t31;
t37 = qJD(2) - t38;
t17 = pkin(6) * t37 + t22;
t35 = -t47 * t43 + t44 * t48;
t30 = -V_base(6) * pkin(1) - t35;
t33 = -t50 * t39 + t52 * V_base(6);
t34 = t39 * t52 + t50 * V_base(6);
t20 = -pkin(2) * t33 - pkin(6) * t34 + t30;
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t10 = -t17 * t49 + t51 * t20;
t25 = t34 * t51 + t37 * t49;
t32 = qJD(3) - t33;
t7 = pkin(3) * t32 - qJ(4) * t25 + t10;
t11 = t51 * t17 + t49 * t20;
t24 = -t34 * t49 + t37 * t51;
t9 = qJ(4) * t24 + t11;
t4 = t46 * t7 + t53 * t9;
t21 = t27 * t52 - t50 * t31;
t16 = -t37 * pkin(2) - t21;
t3 = -t46 * t9 + t53 * t7;
t12 = -t24 * pkin(3) + qJD(4) + t16;
t14 = t46 * t24 + t53 * t25;
t13 = -t53 * t24 + t25 * t46;
t5 = t13 * pkin(4) - t14 * qJ(5) + t12;
t2 = qJ(5) * t32 + t4;
t1 = -t32 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t35 ^ 2 + t36 ^ 2 + t45 ^ 2) / 0.2e1 + m(3) * (t21 ^ 2 + t22 ^ 2 + t30 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t11 ^ 2 + t16 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t45 * mrSges(2,2) - t35 * mrSges(2,3) + Ifges(2,1) * t39 / 0.2e1) * t39 + (t21 * mrSges(3,1) - t22 * mrSges(3,2) + Ifges(3,3) * t37 / 0.2e1) * t37 + (t16 * mrSges(4,2) - t10 * mrSges(4,3) + Ifges(4,1) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t45 * mrSges(2,1) + t36 * mrSges(2,3) + Ifges(2,4) * t39 + Ifges(2,2) * t38 / 0.2e1) * t38 + (t30 * mrSges(3,2) - t21 * mrSges(3,3) + Ifges(3,5) * t37 + Ifges(3,1) * t34 / 0.2e1) * t34 + (-t16 * mrSges(4,1) + t11 * mrSges(4,3) + Ifges(4,4) * t25 + Ifges(4,2) * t24 / 0.2e1) * t24 + (-t30 * mrSges(3,1) + t22 * mrSges(3,3) + Ifges(3,4) * t34 + Ifges(3,6) * t37 + Ifges(3,2) * t33 / 0.2e1) * t33 + (t12 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t14) * t14 + (t12 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t13 + (-Ifges(5,4) + Ifges(6,5)) * t14) * t13 + (V_base(2) * mrSges(1,1) + t35 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t36 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t39 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t38 + (Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6)) * V_base(6) + (t10 * mrSges(4,1) + t3 * mrSges(5,1) - t1 * mrSges(6,1) - t11 * mrSges(4,2) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + Ifges(4,5) * t25 + Ifges(4,6) * t24 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t32 + (Ifges(6,4) + Ifges(5,5)) * t14 + (-Ifges(5,6) + Ifges(6,6)) * t13) * t32;
T = t6;

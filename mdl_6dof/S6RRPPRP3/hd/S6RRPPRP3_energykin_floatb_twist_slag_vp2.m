% Calculate kinetic energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:28
% EndTime: 2019-03-09 08:33:29
% DurationCPUTime: 0.70s
% Computational Cost: add. (1275->147), mult. (1581->179), div. (0->0), fcn. (1092->6), ass. (0->44)
t60 = -pkin(3) - pkin(8);
t52 = sin(qJ(1));
t59 = cos(qJ(1));
t38 = -t52 * V_base(4) + t59 * V_base(5);
t37 = qJD(2) - t38;
t39 = t52 * V_base(5) + t59 * V_base(4);
t49 = V_base(6) + qJD(1);
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t33 = t39 * t54 + t49 * t51;
t23 = -pkin(1) * t38 - pkin(7) * t39 + V_base(3);
t45 = V_base(5) * pkin(6) + V_base(1);
t46 = -V_base(4) * pkin(6) + V_base(2);
t35 = t59 * t45 + t52 * t46;
t27 = pkin(7) * t49 + t35;
t18 = t23 * t54 - t51 * t27;
t58 = qJD(3) - t18;
t56 = -qJ(4) * t33 + t58;
t10 = (-pkin(2) + t60) * t37 + t56;
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t32 = t39 * t51 - t54 * t49;
t34 = -t52 * t45 + t59 * t46;
t26 = -t49 * pkin(1) - t34;
t17 = t32 * pkin(2) - t33 * qJ(3) + t26;
t57 = qJD(4) - t17;
t8 = pkin(4) * t33 + t32 * t60 + t57;
t4 = t53 * t10 + t50 * t8;
t19 = t51 * t23 + t54 * t27;
t16 = t37 * qJ(3) + t19;
t3 = -t10 * t50 + t53 * t8;
t13 = -t32 * qJ(4) - t16;
t12 = pkin(4) * t37 - t13;
t55 = V_base(3) ^ 2;
t31 = qJD(5) + t33;
t21 = t32 * t53 - t37 * t50;
t20 = -t32 * t50 - t37 * t53;
t15 = -pkin(2) * t37 + t58;
t14 = -pkin(3) * t32 + t57;
t11 = (-pkin(2) - pkin(3)) * t37 + t56;
t5 = -pkin(5) * t20 + qJD(6) + t12;
t2 = qJ(6) * t20 + t4;
t1 = pkin(5) * t31 - qJ(6) * t21 + t3;
t6 = m(2) * (t34 ^ 2 + t35 ^ 2 + t55) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t55) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t26 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t13 ^ 2 + t14 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t17 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t34 * mrSges(2,1) - t35 * mrSges(2,2) + Ifges(2,3) * t49 / 0.2e1) * t49 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t34 * mrSges(2,3) + Ifges(2,5) * t49 + Ifges(2,1) * t39 / 0.2e1) * t39 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t35 * mrSges(2,3) + Ifges(2,4) * t39 + Ifges(2,6) * t49 + Ifges(2,2) * t38 / 0.2e1) * t38 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t31) * t31 + (t12 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t21 + (Ifges(6,5) + Ifges(7,5)) * t31) * t21 + (-t12 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t20 + (Ifges(6,6) + Ifges(7,6)) * t31 + (Ifges(6,4) + Ifges(7,4)) * t21) * t20 + (t18 * mrSges(3,1) - t15 * mrSges(4,1) - t13 * mrSges(5,1) - t19 * mrSges(3,2) + t11 * mrSges(5,2) + t16 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t37) * t37 + (t14 * mrSges(5,1) + t26 * mrSges(3,2) + t15 * mrSges(4,2) - t18 * mrSges(3,3) - t17 * mrSges(4,3) - t11 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t33 + (Ifges(4,4) + Ifges(3,5) + Ifges(5,6)) * t37) * t33 + (t26 * mrSges(3,1) + t17 * mrSges(4,1) - t16 * mrSges(4,2) + t14 * mrSges(5,2) - t19 * mrSges(3,3) - t13 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t32 + (-Ifges(5,5) - Ifges(3,6) + Ifges(4,6)) * t37 + (-Ifges(3,4) - Ifges(5,4) + Ifges(4,5)) * t33) * t32;
T  = t6;

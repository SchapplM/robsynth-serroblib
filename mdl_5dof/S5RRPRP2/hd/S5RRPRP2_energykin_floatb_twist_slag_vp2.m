% Calculate kinetic energy for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:23
% EndTime: 2019-12-31 19:49:23
% DurationCPUTime: 0.73s
% Computational Cost: add. (1605->128), mult. (2392->176), div. (0->0), fcn. (1852->8), ass. (0->43)
t50 = sin(qJ(1));
t52 = cos(qJ(1));
t36 = -t50 * V_base(4) + t52 * V_base(5);
t37 = t50 * V_base(5) + t52 * V_base(4);
t49 = sin(qJ(2));
t51 = cos(qJ(2));
t31 = t36 * t51 - t37 * t49;
t32 = t36 * t49 + t37 * t51;
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t23 = t31 * t47 - t32 * t46;
t24 = t31 * t46 + t32 * t47;
t35 = -pkin(1) * t36 + V_base(3);
t25 = -pkin(2) * t31 + qJD(3) + t35;
t12 = -pkin(3) * t23 - pkin(7) * t24 + t25;
t48 = sin(qJ(4));
t54 = cos(qJ(4));
t42 = V_base(5) * pkin(5) + V_base(1);
t43 = -V_base(4) * pkin(5) + V_base(2);
t33 = -t42 * t50 + t52 * t43;
t45 = V_base(6) + qJD(1);
t28 = pkin(1) * t45 - pkin(6) * t37 + t33;
t34 = t52 * t42 + t50 * t43;
t30 = pkin(6) * t36 + t34;
t20 = t51 * t28 - t30 * t49;
t44 = qJD(2) + t45;
t14 = pkin(2) * t44 - qJ(3) * t32 + t20;
t21 = t49 * t28 + t51 * t30;
t17 = qJ(3) * t31 + t21;
t10 = t46 * t14 + t47 * t17;
t8 = pkin(7) * t44 + t10;
t4 = t48 * t12 + t54 * t8;
t9 = t14 * t47 - t46 * t17;
t7 = -pkin(3) * t44 - t9;
t3 = t54 * t12 - t48 * t8;
t53 = V_base(3) ^ 2;
t22 = qJD(4) - t23;
t19 = t54 * t24 + t48 * t44;
t18 = t24 * t48 - t54 * t44;
t5 = pkin(4) * t18 - qJ(5) * t19 + t7;
t2 = qJ(5) * t22 + t4;
t1 = -t22 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t33 ^ 2 + t34 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t35 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t25 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t45 / 0.2e1) * t45 + (t35 * mrSges(3,2) - t20 * mrSges(3,3) + Ifges(3,1) * t32 / 0.2e1) * t32 + (t25 * mrSges(4,2) - t9 * mrSges(4,3) + Ifges(4,1) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) + Ifges(2,5) * t45 + Ifges(2,1) * t37 / 0.2e1) * t37 + (-t35 * mrSges(3,1) + t21 * mrSges(3,3) + Ifges(3,4) * t32 + Ifges(3,2) * t31 / 0.2e1) * t31 + (-t25 * mrSges(4,1) + t10 * mrSges(4,3) + Ifges(4,4) * t24 + Ifges(4,2) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t37 + Ifges(2,6) * t45 + Ifges(2,2) * t36 / 0.2e1) * t36 + (t3 * mrSges(5,1) - t1 * mrSges(6,1) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t22) * t22 + (t7 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t19 + (Ifges(6,4) + Ifges(5,5)) * t22) * t19 + (t20 * mrSges(3,1) + t9 * mrSges(4,1) - t21 * mrSges(3,2) - t10 * mrSges(4,2) + Ifges(3,5) * t32 + Ifges(4,5) * t24 + Ifges(3,6) * t31 + Ifges(4,6) * t23 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t44) * t44 + (t7 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t18 + (-Ifges(5,6) + Ifges(6,6)) * t22 + (-Ifges(5,4) + Ifges(6,5)) * t19) * t18;
T = t6;

% Calculate kinetic energy for
% S5RRPPR5
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:48
% EndTime: 2019-12-31 19:28:49
% DurationCPUTime: 0.86s
% Computational Cost: add. (1399->129), mult. (1846->176), div. (0->0), fcn. (1372->8), ass. (0->46)
t60 = -pkin(3) - pkin(4);
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t37 = -t50 * V_base(4) + t53 * V_base(5);
t38 = t50 * V_base(5) + t53 * V_base(4);
t25 = -pkin(1) * t37 - pkin(6) * t38 + V_base(3);
t43 = V_base(5) * pkin(5) + V_base(1);
t44 = -V_base(4) * pkin(5) + V_base(2);
t33 = t53 * t43 + t50 * t44;
t46 = V_base(6) + qJD(1);
t29 = pkin(6) * t46 + t33;
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t18 = t52 * t25 - t29 * t49;
t31 = t38 * t52 + t46 * t49;
t36 = qJD(2) - t37;
t14 = pkin(2) * t36 - qJ(3) * t31 + t18;
t19 = t49 * t25 + t52 * t29;
t30 = -t38 * t49 + t46 * t52;
t17 = qJ(3) * t30 + t19;
t47 = sin(pkin(8));
t59 = cos(pkin(8));
t9 = t47 * t14 + t59 * t17;
t32 = -t50 * t43 + t53 * t44;
t7 = t36 * qJ(4) + t9;
t58 = pkin(1) * t46 + t32;
t8 = t59 * t14 - t47 * t17;
t57 = qJD(4) - t8;
t56 = pkin(2) * t30 - qJD(3) + t58;
t21 = t47 * t30 + t59 * t31;
t55 = qJ(4) * t21 + t56;
t54 = V_base(3) ^ 2;
t51 = cos(qJ(5));
t48 = sin(qJ(5));
t35 = qJD(5) - t36;
t20 = -t59 * t30 + t31 * t47;
t12 = t20 * t48 + t21 * t51;
t11 = t20 * t51 - t21 * t48;
t10 = pkin(3) * t20 - t55;
t6 = -t36 * pkin(3) + t57;
t5 = t60 * t20 + t55;
t4 = pkin(7) * t20 + t7;
t3 = -t21 * pkin(7) + t60 * t36 + t57;
t2 = t3 * t48 + t4 * t51;
t1 = t3 * t51 - t4 * t48;
t13 = m(2) * (t32 ^ 2 + t33 ^ 2 + t54) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t54) / 0.2e1 + m(3) * (t18 ^ 2 + t19 ^ 2 + t58 ^ 2) / 0.2e1 + m(4) * (t56 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t32 * mrSges(2,1) - t33 * mrSges(2,2) + Ifges(2,3) * t46 / 0.2e1) * t46 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t35 / 0.2e1) * t35 + (-t58 * mrSges(3,2) - t18 * mrSges(3,3) + Ifges(3,1) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t32 * mrSges(2,3) + Ifges(2,5) * t46 + Ifges(2,1) * t38 / 0.2e1) * t38 + (t58 * mrSges(3,1) + t19 * mrSges(3,3) + Ifges(3,4) * t31 + Ifges(3,2) * t30 / 0.2e1) * t30 + (t5 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t35 + Ifges(6,1) * t12 / 0.2e1) * t12 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t33 * mrSges(2,3) + Ifges(2,4) * t38 + Ifges(2,6) * t46 + Ifges(2,2) * t37 / 0.2e1) * t37 + (-t5 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t12 + Ifges(6,6) * t35 + Ifges(6,2) * t11 / 0.2e1) * t11 + (-t56 * mrSges(4,2) + t6 * mrSges(5,2) - t8 * mrSges(4,3) - t10 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t21) * t21 + (-t56 * mrSges(4,1) + t10 * mrSges(5,1) - t7 * mrSges(5,2) - t9 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t20 + (-Ifges(4,4) + Ifges(5,5)) * t21) * t20 + (t18 * mrSges(3,1) + t8 * mrSges(4,1) - t6 * mrSges(5,1) - t19 * mrSges(3,2) - t9 * mrSges(4,2) + t7 * mrSges(5,3) + Ifges(3,5) * t31 + Ifges(3,6) * t30 + (Ifges(3,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t36 + (Ifges(5,4) + Ifges(4,5)) * t21 + (-Ifges(4,6) + Ifges(5,6)) * t20) * t36;
T = t13;

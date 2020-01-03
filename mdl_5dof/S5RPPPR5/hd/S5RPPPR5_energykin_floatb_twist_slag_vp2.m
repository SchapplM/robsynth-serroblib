% Calculate kinetic energy for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:18
% EndTime: 2019-12-31 17:46:19
% DurationCPUTime: 0.66s
% Computational Cost: add. (1143->129), mult. (1566->174), div. (0->0), fcn. (1092->8), ass. (0->43)
t51 = sin(qJ(1));
t56 = cos(qJ(1));
t37 = t51 * V_base(4) - t56 * V_base(5);
t38 = t51 * V_base(5) + t56 * V_base(4);
t29 = t37 * pkin(1) - t38 * qJ(2) + V_base(3);
t20 = -pkin(2) * t37 + qJD(3) - t29;
t48 = sin(pkin(7));
t55 = cos(pkin(7));
t27 = -t55 * t37 + t38 * t48;
t28 = t48 * t37 + t55 * t38;
t10 = pkin(3) * t27 - qJ(4) * t28 + t20;
t46 = V_base(6) + qJD(1);
t42 = V_base(5) * pkin(5) + V_base(1);
t43 = -V_base(4) * pkin(5) + V_base(2);
t32 = -t51 * t42 + t56 * t43;
t54 = qJD(2) - t32;
t19 = -t38 * qJ(3) + (-pkin(1) - pkin(2)) * t46 + t54;
t33 = t56 * t42 + t51 * t43;
t31 = t46 * qJ(2) + t33;
t25 = qJ(3) * t37 + t31;
t15 = t48 * t19 + t55 * t25;
t13 = -qJ(4) * t46 + t15;
t47 = sin(pkin(8));
t49 = cos(pkin(8));
t6 = t47 * t10 + t49 * t13;
t5 = t49 * t10 - t13 * t47;
t14 = t55 * t19 - t48 * t25;
t12 = t46 * pkin(3) + qJD(4) - t14;
t53 = V_base(3) ^ 2;
t52 = cos(qJ(5));
t50 = sin(qJ(5));
t30 = -t46 * pkin(1) + t54;
t26 = qJD(5) + t27;
t22 = t28 * t49 - t46 * t47;
t21 = -t28 * t47 - t46 * t49;
t17 = t21 * t50 + t22 * t52;
t16 = t21 * t52 - t22 * t50;
t7 = -t21 * pkin(4) + t12;
t4 = pkin(6) * t21 + t6;
t3 = pkin(4) * t27 - pkin(6) * t22 + t5;
t2 = t3 * t50 + t4 * t52;
t1 = t3 * t52 - t4 * t50;
t8 = m(2) * (t32 ^ 2 + t33 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t20 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t20 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,1) * t28 / 0.2e1) * t28 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t26 / 0.2e1) * t26 + (t12 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t22 / 0.2e1) * t22 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t12 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t22 + Ifges(5,2) * t21 / 0.2e1) * t21 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t26 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t26 + Ifges(6,2) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(2,2) + t30 * mrSges(3,2) - t32 * mrSges(2,3) - t29 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t38) * t38 + (V_base(3) * mrSges(2,1) + t29 * mrSges(3,1) - t31 * mrSges(3,2) - t33 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t37 + (-Ifges(2,4) + Ifges(3,5)) * t38) * t37 + (t20 * mrSges(4,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) - t15 * mrSges(4,3) - Ifges(4,4) * t28 + Ifges(5,5) * t22 + Ifges(5,6) * t21 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t27) * t27 + (t32 * mrSges(2,1) - t30 * mrSges(3,1) - t14 * mrSges(4,1) - t33 * mrSges(2,2) + t15 * mrSges(4,2) + t31 * mrSges(3,3) - Ifges(4,5) * t28 + Ifges(4,6) * t27 + (Ifges(4,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t46 + (Ifges(3,4) + Ifges(2,5)) * t38 + (-Ifges(2,6) + Ifges(3,6)) * t37) * t46;
T = t8;

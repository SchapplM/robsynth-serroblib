% Calculate kinetic energy for
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:31
% EndTime: 2019-12-31 20:57:32
% DurationCPUTime: 0.73s
% Computational Cost: add. (1137->125), mult. (1468->161), div. (0->0), fcn. (1056->6), ass. (0->40)
t54 = -pkin(3) - pkin(4);
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t33 = -t45 * V_base(4) + t47 * V_base(5);
t34 = t45 * V_base(5) + t47 * V_base(4);
t21 = -pkin(1) * t33 - pkin(6) * t34 + V_base(3);
t39 = V_base(5) * pkin(5) + V_base(1);
t40 = -V_base(4) * pkin(5) + V_base(2);
t29 = t47 * t39 + t45 * t40;
t42 = V_base(6) + qJD(1);
t25 = pkin(6) * t42 + t29;
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t14 = t46 * t21 - t25 * t44;
t27 = t34 * t46 + t42 * t44;
t32 = qJD(2) - t33;
t10 = pkin(2) * t32 - pkin(7) * t27 + t14;
t15 = t44 * t21 + t46 * t25;
t26 = -t34 * t44 + t42 * t46;
t13 = pkin(7) * t26 + t15;
t43 = sin(qJ(3));
t53 = cos(qJ(3));
t7 = t43 * t10 + t53 * t13;
t28 = -t45 * t39 + t47 * t40;
t31 = qJD(3) + t32;
t5 = t31 * qJ(4) + t7;
t52 = pkin(1) * t42 + t28;
t6 = t53 * t10 - t43 * t13;
t51 = qJD(4) - t6;
t50 = pkin(2) * t26 + t52;
t17 = t43 * t26 + t53 * t27;
t49 = qJ(4) * t17 + t50;
t48 = V_base(3) ^ 2;
t16 = -t53 * t26 + t27 * t43;
t8 = pkin(3) * t16 - t49;
t4 = -t31 * pkin(3) + t51;
t3 = t54 * t16 + qJD(5) + t49;
t2 = qJ(5) * t16 + t5;
t1 = -t17 * qJ(5) + t54 * t31 + t51;
t9 = m(2) * (t28 ^ 2 + t29 ^ 2 + t48) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t48) / 0.2e1 + m(4) * (t50 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(3) * (t14 ^ 2 + t15 ^ 2 + t52 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t28 * mrSges(2,1) - t29 * mrSges(2,2) + Ifges(2,3) * t42 / 0.2e1) * t42 + (t14 * mrSges(3,1) - t15 * mrSges(3,2) + Ifges(3,3) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t28 * mrSges(2,3) + Ifges(2,5) * t42 + Ifges(2,1) * t34 / 0.2e1) * t34 + (-t52 * mrSges(3,2) - t14 * mrSges(3,3) + Ifges(3,5) * t32 + Ifges(3,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t29 * mrSges(2,3) + Ifges(2,4) * t34 + Ifges(2,6) * t42 + Ifges(2,2) * t33 / 0.2e1) * t33 + (t52 * mrSges(3,1) + t15 * mrSges(3,3) + Ifges(3,4) * t27 + Ifges(3,6) * t32 + Ifges(3,2) * t26 / 0.2e1) * t26 + (t6 * mrSges(4,1) - t4 * mrSges(5,1) - t1 * mrSges(6,1) - t7 * mrSges(4,2) + t2 * mrSges(6,2) + t5 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t31) * t31 + (-t50 * mrSges(4,2) + t4 * mrSges(5,2) + t3 * mrSges(6,2) - t6 * mrSges(4,3) - t8 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t17 + (Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t31) * t17 + (-t50 * mrSges(4,1) + t8 * mrSges(5,1) - t3 * mrSges(6,1) - t5 * mrSges(5,2) - t7 * mrSges(4,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t16 + (-Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t31 + (-Ifges(4,4) + Ifges(6,4) + Ifges(5,5)) * t17) * t16;
T = t9;

% Calculate kinetic energy for
% S5RPRRR13
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR13_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:27
% EndTime: 2019-12-31 19:14:28
% DurationCPUTime: 0.77s
% Computational Cost: add. (1203->129), mult. (1512->178), div. (0->0), fcn. (1044->8), ass. (0->47)
t58 = pkin(1) + pkin(6);
t50 = sin(qJ(1));
t57 = cos(qJ(1));
t37 = t50 * V_base(5) + t57 * V_base(4);
t46 = V_base(6) + qJD(1);
t42 = V_base(5) * pkin(5) + V_base(1);
t43 = -V_base(4) * pkin(5) + V_base(2);
t33 = -t50 * t42 + t57 * t43;
t55 = qJD(2) - t33;
t19 = t37 * pkin(2) - t58 * t46 + t55;
t36 = t50 * V_base(4) - t57 * V_base(5);
t56 = -qJ(2) * t37 + V_base(3);
t24 = t58 * t36 + t56;
t49 = sin(qJ(3));
t53 = cos(qJ(3));
t12 = t49 * t19 + t53 * t24;
t35 = qJD(3) + t37;
t10 = pkin(7) * t35 + t12;
t34 = t57 * t42 + t50 * t43;
t28 = -t46 * qJ(2) - t34;
t25 = -pkin(2) * t36 - t28;
t31 = t36 * t53 - t49 * t46;
t32 = t36 * t49 + t46 * t53;
t17 = -pkin(3) * t31 - pkin(7) * t32 + t25;
t48 = sin(qJ(4));
t52 = cos(qJ(4));
t6 = t52 * t10 + t48 * t17;
t5 = -t10 * t48 + t52 * t17;
t11 = t19 * t53 - t49 * t24;
t30 = qJD(4) - t31;
t9 = -pkin(3) * t35 - t11;
t54 = V_base(3) ^ 2;
t51 = cos(qJ(5));
t47 = sin(qJ(5));
t29 = qJD(5) + t30;
t27 = -t46 * pkin(1) + t55;
t26 = pkin(1) * t36 + t56;
t23 = t32 * t52 + t35 * t48;
t22 = -t32 * t48 + t35 * t52;
t14 = t22 * t47 + t23 * t51;
t13 = t22 * t51 - t23 * t47;
t7 = -pkin(4) * t22 + t9;
t4 = pkin(8) * t22 + t6;
t3 = pkin(4) * t30 - pkin(8) * t23 + t5;
t2 = t3 * t47 + t4 * t51;
t1 = t3 * t51 - t4 * t47;
t8 = m(2) * (t33 ^ 2 + t34 ^ 2 + t54) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t54) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t25 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t35 / 0.2e1) * t35 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t30 / 0.2e1) * t30 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t25 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t35 + Ifges(4,1) * t32 / 0.2e1) * t32 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t30 + Ifges(5,1) * t23 / 0.2e1) * t23 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t25 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t35 + Ifges(4,2) * t31 / 0.2e1) * t31 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t23 + Ifges(5,6) * t30 + Ifges(5,2) * t22 / 0.2e1) * t22 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t29 + Ifges(6,2) * t13 / 0.2e1) * t13 + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + t27 * mrSges(3,2) - t28 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t46) * t46 + (t27 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) - t26 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t37 + (-Ifges(3,4) + Ifges(2,5)) * t46) * t37 + (V_base(3) * mrSges(2,1) + t28 * mrSges(3,1) - t26 * mrSges(3,2) - t34 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t36 + (Ifges(3,5) - Ifges(2,6)) * t46 + (-Ifges(2,4) - Ifges(3,6)) * t37) * t36;
T = t8;

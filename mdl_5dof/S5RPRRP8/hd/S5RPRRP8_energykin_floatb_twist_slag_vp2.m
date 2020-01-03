% Calculate kinetic energy for
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:04
% EndTime: 2019-12-31 18:47:05
% DurationCPUTime: 0.65s
% Computational Cost: add. (931->125), mult. (1204->161), div. (0->0), fcn. (796->6), ass. (0->38)
t45 = sin(qJ(1));
t50 = cos(qJ(1));
t32 = t45 * V_base(5) + t50 * V_base(4);
t42 = V_base(6) + qJD(1);
t37 = V_base(5) * pkin(5) + V_base(1);
t38 = -V_base(4) * pkin(5) + V_base(2);
t27 = -t45 * t37 + t50 * t38;
t48 = qJD(2) - t27;
t14 = -t32 * pkin(6) + (-pkin(1) - pkin(2)) * t42 + t48;
t28 = t50 * t37 + t45 * t38;
t26 = t42 * qJ(2) + t28;
t31 = t45 * V_base(4) - t50 * V_base(5);
t20 = pkin(6) * t31 + t26;
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t12 = t44 * t14 + t46 * t20;
t40 = qJD(3) - t42;
t10 = pkin(7) * t40 + t12;
t43 = sin(qJ(4));
t49 = cos(qJ(4));
t24 = t31 * pkin(1) - t32 * qJ(2) + V_base(3);
t17 = -pkin(2) * t31 - t24;
t22 = t31 * t46 - t32 * t44;
t23 = t31 * t44 + t32 * t46;
t7 = -pkin(3) * t22 - pkin(7) * t23 + t17;
t4 = t49 * t10 + t43 * t7;
t11 = t14 * t46 - t44 * t20;
t9 = -pkin(3) * t40 - t11;
t3 = -t43 * t10 + t49 * t7;
t47 = V_base(3) ^ 2;
t25 = -t42 * pkin(1) + t48;
t21 = qJD(4) - t22;
t16 = t49 * t23 + t43 * t40;
t15 = t23 * t43 - t49 * t40;
t5 = pkin(4) * t15 - qJ(5) * t16 + t9;
t2 = qJ(5) * t21 + t4;
t1 = -t21 * pkin(4) + qJD(5) - t3;
t6 = m(2) * (t27 ^ 2 + t28 ^ 2 + t47) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t47) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t17 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t40 / 0.2e1) * t40 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t17 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t40 + Ifges(4,1) * t23 / 0.2e1) * t23 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t17 * mrSges(4,1) + t12 * mrSges(4,3) + Ifges(4,4) * t23 + Ifges(4,6) * t40 + Ifges(4,2) * t22 / 0.2e1) * t22 + (t27 * mrSges(2,1) - t25 * mrSges(3,1) - t28 * mrSges(2,2) + t26 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t42) * t42 + (t3 * mrSges(5,1) - t1 * mrSges(6,1) - t4 * mrSges(5,2) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t21) * t21 + (V_base(3) * mrSges(2,2) + t25 * mrSges(3,2) - t27 * mrSges(2,3) - t24 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t32 + (Ifges(3,4) + Ifges(2,5)) * t42) * t32 + (t9 * mrSges(5,2) + t1 * mrSges(6,2) - t3 * mrSges(5,3) - t5 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t16 + (Ifges(6,4) + Ifges(5,5)) * t21) * t16 + (V_base(3) * mrSges(2,1) + t24 * mrSges(3,1) - t26 * mrSges(3,2) - t28 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t31 + (-Ifges(2,6) + Ifges(3,6)) * t42 + (-Ifges(2,4) + Ifges(3,5)) * t32) * t31 + (t9 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(6,2) - t4 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t15 + (-Ifges(5,6) + Ifges(6,6)) * t21 + (-Ifges(5,4) + Ifges(6,5)) * t16) * t15;
T = t6;

% Calculate kinetic energy for
% S5RPRRP12
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:11
% EndTime: 2019-12-31 18:56:12
% DurationCPUTime: 0.65s
% Computational Cost: add. (929->125), mult. (1166->161), div. (0->0), fcn. (764->6), ass. (0->40)
t50 = pkin(1) + pkin(6);
t36 = V_base(5) * pkin(5) + V_base(1);
t37 = -V_base(4) * pkin(5) + V_base(2);
t43 = sin(qJ(1));
t49 = cos(qJ(1));
t29 = t49 * t36 + t43 * t37;
t40 = V_base(6) + qJD(1);
t24 = -t40 * qJ(2) - t29;
t31 = t43 * V_base(4) - t49 * V_base(5);
t21 = -pkin(2) * t31 - t24;
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t26 = t31 * t45 - t40 * t42;
t27 = t31 * t42 + t40 * t45;
t13 = -pkin(3) * t26 - pkin(7) * t27 + t21;
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t32 = t43 * V_base(5) + t49 * V_base(4);
t28 = -t43 * t36 + t49 * t37;
t47 = qJD(2) - t28;
t15 = t32 * pkin(2) - t50 * t40 + t47;
t48 = -qJ(2) * t32 + V_base(3);
t20 = t50 * t31 + t48;
t10 = t42 * t15 + t45 * t20;
t30 = qJD(3) + t32;
t8 = pkin(7) * t30 + t10;
t4 = t41 * t13 + t44 * t8;
t3 = t44 * t13 - t41 * t8;
t9 = t15 * t45 - t42 * t20;
t7 = -pkin(3) * t30 - t9;
t46 = V_base(3) ^ 2;
t25 = qJD(4) - t26;
t23 = -t40 * pkin(1) + t47;
t22 = pkin(1) * t31 + t48;
t19 = t27 * t44 + t30 * t41;
t18 = -t27 * t41 + t30 * t44;
t5 = -pkin(4) * t18 + qJD(5) + t7;
t2 = qJ(5) * t18 + t4;
t1 = pkin(4) * t25 - qJ(5) * t19 + t3;
t6 = m(2) * (t28 ^ 2 + t29 ^ 2 + t46) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t46) / 0.2e1 + m(4) * (t10 ^ 2 + t21 ^ 2 + t9 ^ 2) / 0.2e1 + m(3) * (t22 ^ 2 + t23 ^ 2 + t24 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t9 * mrSges(4,1) - t10 * mrSges(4,2) + Ifges(4,3) * t30 / 0.2e1) * t30 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t21 * mrSges(4,2) - t9 * mrSges(4,3) + Ifges(4,5) * t30 + Ifges(4,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t21 * mrSges(4,1) + t10 * mrSges(4,3) + Ifges(4,4) * t27 + Ifges(4,6) * t30 + Ifges(4,2) * t26 / 0.2e1) * t26 + (t28 * mrSges(2,1) - t29 * mrSges(2,2) + t23 * mrSges(3,2) - t24 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t40) * t40 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t25) * t25 + (t23 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t28 * mrSges(2,3) - t22 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t32 + (-Ifges(3,4) + Ifges(2,5)) * t40) * t32 + (t7 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t19 + (Ifges(5,5) + Ifges(6,5)) * t25) * t19 + (V_base(3) * mrSges(2,1) + t24 * mrSges(3,1) - t22 * mrSges(3,2) - t29 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t31 + (Ifges(3,5) - Ifges(2,6)) * t40 + (-Ifges(2,4) - Ifges(3,6)) * t32) * t31 + (-t7 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t18 + (Ifges(5,6) + Ifges(6,6)) * t25 + (Ifges(5,4) + Ifges(6,4)) * t19) * t18;
T = t6;

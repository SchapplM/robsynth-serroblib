% Calculate kinetic energy for
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP4_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:23
% EndTime: 2019-12-31 17:34:24
% DurationCPUTime: 0.59s
% Computational Cost: add. (837->125), mult. (1208->160), div. (0->0), fcn. (800->6), ass. (0->37)
t43 = sin(pkin(7));
t49 = cos(pkin(7));
t33 = t43 * V_base(5) + t49 * V_base(4);
t37 = V_base(5) * qJ(1) + V_base(1);
t38 = -V_base(4) * qJ(1) + V_base(2);
t28 = -t43 * t37 + t38 * t49;
t48 = qJD(2) - t28;
t16 = -t33 * pkin(5) + (-pkin(1) - pkin(2)) * V_base(6) + t48;
t29 = t49 * t37 + t43 * t38;
t27 = V_base(6) * qJ(2) + t29;
t32 = t43 * V_base(4) - t49 * V_base(5);
t21 = pkin(5) * t32 + t27;
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t13 = t45 * t16 + t47 * t21;
t40 = -V_base(6) + qJD(3);
t11 = pkin(6) * t40 + t13;
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t42 = V_base(3) + qJD(1);
t23 = t32 * pkin(1) - t33 * qJ(2) + t42;
t15 = -pkin(2) * t32 - t23;
t24 = t32 * t47 - t45 * t33;
t25 = t45 * t32 + t33 * t47;
t8 = -pkin(3) * t24 - pkin(6) * t25 + t15;
t4 = t46 * t11 + t44 * t8;
t3 = -t11 * t44 + t46 * t8;
t12 = t16 * t47 - t45 * t21;
t10 = -t40 * pkin(3) - t12;
t26 = -V_base(6) * pkin(1) + t48;
t22 = qJD(4) - t24;
t18 = t25 * t46 + t40 * t44;
t17 = -t25 * t44 + t40 * t46;
t5 = -t17 * pkin(4) + qJD(5) + t10;
t2 = qJ(5) * t17 + t4;
t1 = pkin(4) * t22 - qJ(5) * t18 + t3;
t6 = m(2) * (t28 ^ 2 + t29 ^ 2 + t42 ^ 2) / 0.2e1 + m(3) * (t23 ^ 2 + t26 ^ 2 + t27 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t15 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t12 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,3) * t40 / 0.2e1) * t40 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t15 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,5) * t40 + Ifges(4,1) * t25 / 0.2e1) * t25 + (-t15 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t25 + Ifges(4,6) * t40 + Ifges(4,2) * t24 / 0.2e1) * t24 + (t42 * mrSges(2,2) + t26 * mrSges(3,2) - t28 * mrSges(2,3) - t23 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t33) * t33 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t22) * t22 + (t42 * mrSges(2,1) + t23 * mrSges(3,1) - t27 * mrSges(3,2) - t29 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t32 + (-Ifges(2,4) + Ifges(3,5)) * t33) * t32 + (t10 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t18 + (Ifges(5,5) + Ifges(6,5)) * t22) * t18 + (-t10 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t17 + (Ifges(5,6) + Ifges(6,6)) * t22 + (Ifges(5,4) + Ifges(6,4)) * t18) * t17 + (V_base(2) * mrSges(1,1) + t28 * mrSges(2,1) - t26 * mrSges(3,1) - V_base(1) * mrSges(1,2) - t29 * mrSges(2,2) + t27 * mrSges(3,3) + Ifges(1,5) * V_base(4) + Ifges(1,6) * V_base(5) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6) + (Ifges(3,4) + Ifges(2,5)) * t33 + (-Ifges(2,6) + Ifges(3,6)) * t32) * V_base(6);
T = t6;

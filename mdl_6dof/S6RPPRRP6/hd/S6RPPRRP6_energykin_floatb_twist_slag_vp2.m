% Calculate kinetic energy for
% S6RPPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP6_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:07
% EndTime: 2019-03-09 02:10:08
% DurationCPUTime: 0.70s
% Computational Cost: add. (1123->147), mult. (1399->179), div. (0->0), fcn. (892->6), ass. (0->43)
t42 = V_base(5) * pkin(6) + V_base(1);
t43 = -V_base(4) * pkin(6) + V_base(2);
t51 = sin(qJ(1));
t57 = cos(qJ(1));
t30 = -t51 * t42 + t57 * t43;
t48 = V_base(6) + qJD(1);
t25 = -t48 * pkin(1) + qJD(2) - t30;
t37 = t51 * V_base(5) + t57 * V_base(4);
t19 = t37 * pkin(2) - t48 * qJ(3) + t25;
t14 = -pkin(3) * t37 - t19;
t50 = sin(qJ(4));
t52 = cos(qJ(4));
t28 = t37 * t52 - t48 * t50;
t29 = t37 * t50 + t48 * t52;
t12 = -pkin(4) * t28 - pkin(8) * t29 + t14;
t49 = sin(qJ(5));
t56 = cos(qJ(5));
t36 = t51 * V_base(4) - t57 * V_base(5);
t31 = t57 * t42 + t51 * t43;
t26 = -t48 * qJ(2) - t31;
t55 = qJD(3) - t26;
t15 = -pkin(7) * t48 + (-pkin(2) - pkin(3)) * t36 + t55;
t33 = t36 * qJ(3);
t54 = pkin(1) * t36 + V_base(3);
t18 = t33 + (pkin(7) - qJ(2)) * t37 + t54;
t10 = t50 * t15 + t52 * t18;
t35 = qJD(4) - t36;
t8 = pkin(8) * t35 + t10;
t4 = t49 * t12 + t56 * t8;
t9 = t15 * t52 - t50 * t18;
t7 = -pkin(4) * t35 - t9;
t3 = t12 * t56 - t49 * t8;
t24 = -qJ(2) * t37 + t54;
t53 = V_base(3) ^ 2;
t27 = qJD(5) - t28;
t23 = -pkin(2) * t36 + t55;
t22 = t24 + t33;
t21 = t29 * t56 + t49 * t35;
t20 = t29 * t49 - t35 * t56;
t5 = pkin(5) * t20 - qJ(6) * t21 + t7;
t2 = qJ(6) * t27 + t4;
t1 = -t27 * pkin(5) + qJD(6) - t3;
t6 = m(2) * (t30 ^ 2 + t31 ^ 2 + t53) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t53) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(4) * (t19 ^ 2 + t22 ^ 2 + t23 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t10 ^ 2 + t14 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t9 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t35 / 0.2e1) * t35 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t14 * mrSges(5,2) - t9 * mrSges(5,3) + Ifges(5,5) * t35 + Ifges(5,1) * t29 / 0.2e1) * t29 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t14 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t29 + Ifges(5,6) * t35 + Ifges(5,2) * t28 / 0.2e1) * t28 + (t3 * mrSges(6,1) - t1 * mrSges(7,1) - t4 * mrSges(6,2) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t27) * t27 + (t7 * mrSges(6,2) + t1 * mrSges(7,2) - t3 * mrSges(6,3) - t5 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t21 + (Ifges(7,4) + Ifges(6,5)) * t27) * t21 + (t7 * mrSges(6,1) + t5 * mrSges(7,1) - t2 * mrSges(7,2) - t4 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t20 + (-Ifges(6,6) + Ifges(7,6)) * t27 + (-Ifges(6,4) + Ifges(7,5)) * t21) * t20 + (t30 * mrSges(2,1) - t31 * mrSges(2,2) + t25 * mrSges(3,2) + t23 * mrSges(4,2) - t26 * mrSges(3,3) - t19 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t48) * t48 + (t25 * mrSges(3,1) + t19 * mrSges(4,1) + V_base(3) * mrSges(2,2) - t22 * mrSges(4,2) - t30 * mrSges(2,3) - t24 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t37 + (-Ifges(3,4) + Ifges(2,5) + Ifges(4,5)) * t48) * t37 + (V_base(3) * mrSges(2,1) + t26 * mrSges(3,1) - t23 * mrSges(4,1) - t24 * mrSges(3,2) - t31 * mrSges(2,3) + t22 * mrSges(4,3) + (Ifges(2,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t36 + (Ifges(4,4) + Ifges(3,5) - Ifges(2,6)) * t48 + (-Ifges(2,4) - Ifges(3,6) + Ifges(4,6)) * t37) * t36;
T  = t6;

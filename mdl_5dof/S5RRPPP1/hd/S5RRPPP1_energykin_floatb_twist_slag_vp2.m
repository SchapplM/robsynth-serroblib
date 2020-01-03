% Calculate kinetic energy for
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:16
% EndTime: 2019-12-31 19:23:17
% DurationCPUTime: 0.71s
% Computational Cost: add. (1713->132), mult. (2366->175), div. (0->0), fcn. (1834->8), ass. (0->47)
t59 = pkin(3) + qJ(5);
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t36 = -t48 * V_base(4) + t50 * V_base(5);
t37 = t48 * V_base(5) + t50 * V_base(4);
t26 = -pkin(1) * t36 - pkin(7) * t37 + V_base(3);
t41 = V_base(5) * pkin(6) + V_base(1);
t42 = -V_base(4) * pkin(6) + V_base(2);
t33 = t50 * t41 + t48 * t42;
t43 = V_base(6) + qJD(1);
t29 = pkin(7) * t43 + t33;
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t21 = t47 * t26 + t49 * t29;
t31 = t37 * t49 + t43 * t47;
t58 = qJ(3) * t31;
t57 = cos(pkin(8));
t30 = -t37 * t47 + t43 * t49;
t35 = qJD(2) - t36;
t45 = sin(pkin(5));
t46 = cos(pkin(5));
t54 = t30 * t46 + t35 * t45;
t13 = t54 * qJ(3) + t21;
t20 = t49 * t26 - t29 * t47;
t14 = pkin(2) * t35 - t46 * t58 + t20;
t32 = -t48 * t41 + t42 * t50;
t28 = -pkin(1) * t43 - t32;
t17 = -pkin(2) * t30 - t45 * t58 + t28;
t44 = sin(pkin(8));
t8 = t57 * t13 + (t14 * t46 + t17 * t45) * t44;
t56 = t45 * t57;
t55 = t46 * t57;
t9 = -t14 * t45 + t46 * t17 + qJD(3);
t22 = -t30 * t45 + t35 * t46;
t5 = -qJ(4) * t22 - t8;
t19 = t57 * t31 + t54 * t44;
t53 = -qJ(4) * t19 + t9;
t7 = -t44 * t13 + t14 * t55 + t17 * t56;
t52 = qJD(4) - t7;
t51 = V_base(3) ^ 2;
t18 = -t30 * t55 + t31 * t44 - t35 * t56;
t6 = pkin(3) * t18 + t53;
t4 = -t22 * pkin(3) + t52;
t3 = t59 * t18 + t53;
t2 = -pkin(4) * t18 + qJD(5) - t5;
t1 = t19 * pkin(4) - t59 * t22 + t52;
t10 = m(2) * (t32 ^ 2 + t33 ^ 2 + t51) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t51) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t28 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t7 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t32 * mrSges(2,1) - t33 * mrSges(2,2) + Ifges(2,3) * t43 / 0.2e1) * t43 + (t20 * mrSges(3,1) - t21 * mrSges(3,2) + Ifges(3,3) * t35 / 0.2e1) * t35 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t32 * mrSges(2,3) + Ifges(2,5) * t43 + Ifges(2,1) * t37 / 0.2e1) * t37 + (t28 * mrSges(3,2) - t20 * mrSges(3,3) + Ifges(3,5) * t35 + Ifges(3,1) * t31 / 0.2e1) * t31 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t33 * mrSges(2,3) + Ifges(2,4) * t37 + Ifges(2,6) * t43 + Ifges(2,2) * t36 / 0.2e1) * t36 + (-t28 * mrSges(3,1) + t21 * mrSges(3,3) + Ifges(3,4) * t31 + Ifges(3,6) * t35 + Ifges(3,2) * t30 / 0.2e1) * t30 + (t7 * mrSges(4,1) - t8 * mrSges(4,2) + t4 * mrSges(5,2) + t2 * mrSges(6,2) - t5 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t22) * t22 + (t4 * mrSges(5,1) + t1 * mrSges(6,1) + t9 * mrSges(4,2) - t3 * mrSges(6,2) - t7 * mrSges(4,3) - t6 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t19 + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t22) * t19 + (t9 * mrSges(4,1) + t5 * mrSges(5,1) - t2 * mrSges(6,1) - t6 * mrSges(5,2) - t8 * mrSges(4,3) + t3 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t18 + (Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t22 + (-Ifges(4,4) - Ifges(5,6) + Ifges(6,6)) * t19) * t18;
T = t10;

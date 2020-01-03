% Calculate kinetic energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:03
% EndTime: 2019-12-31 22:08:04
% DurationCPUTime: 0.96s
% Computational Cost: add. (2473->133), mult. (3706->188), div. (0->0), fcn. (3022->10), ass. (0->50)
t47 = pkin(6) * V_base(5) + V_base(1);
t48 = -pkin(6) * V_base(4) + V_base(2);
t55 = sin(qJ(1));
t59 = cos(qJ(1));
t40 = -t47 * t55 + t59 * t48;
t49 = V_base(6) + qJD(1);
t51 = cos(pkin(5));
t43 = t55 * V_base(5) + t59 * V_base(4);
t66 = pkin(7) * t43;
t35 = pkin(1) * t49 - t51 * t66 + t40;
t42 = -t55 * V_base(4) + t59 * V_base(5);
t50 = sin(pkin(5));
t38 = -pkin(1) * t42 - t50 * t66 + V_base(3);
t67 = t35 * t51 + t38 * t50;
t41 = t59 * t47 + t55 * t48;
t61 = t42 * t51 + t49 * t50;
t32 = pkin(7) * t61 + t41;
t54 = sin(qJ(2));
t58 = cos(qJ(2));
t20 = -t54 * t32 + t58 * t67;
t33 = -t43 * t54 + t58 * t61;
t39 = -t42 * t50 + t49 * t51 + qJD(2);
t18 = -pkin(2) * t39 - t20;
t34 = t43 * t58 + t54 * t61;
t53 = sin(qJ(3));
t57 = cos(qJ(3));
t26 = -t34 * t53 + t39 * t57;
t27 = t34 * t57 + t39 * t53;
t13 = -pkin(3) * t26 - pkin(9) * t27 + t18;
t52 = sin(qJ(4));
t56 = cos(qJ(4));
t24 = -t35 * t50 + t51 * t38;
t15 = -pkin(2) * t33 - pkin(8) * t34 + t24;
t21 = t58 * t32 + t54 * t67;
t19 = pkin(8) * t39 + t21;
t10 = t53 * t15 + t57 * t19;
t31 = qJD(3) - t33;
t8 = pkin(9) * t31 + t10;
t4 = t52 * t13 + t56 * t8;
t3 = t56 * t13 - t52 * t8;
t9 = t15 * t57 - t53 * t19;
t7 = -pkin(3) * t31 - t9;
t60 = V_base(3) ^ 2;
t25 = qJD(4) - t26;
t23 = t27 * t56 + t31 * t52;
t22 = -t27 * t52 + t31 * t56;
t5 = -pkin(4) * t22 + qJD(5) + t7;
t2 = qJ(5) * t22 + t4;
t1 = pkin(4) * t25 - qJ(5) * t23 + t3;
t6 = m(2) * (t40 ^ 2 + t41 ^ 2 + t60) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t60) / 0.2e1 + m(3) * (t20 ^ 2 + t21 ^ 2 + t24 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t18 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t40 * mrSges(2,1) - t41 * mrSges(2,2) + Ifges(2,3) * t49 / 0.2e1) * t49 + (t20 * mrSges(3,1) - t21 * mrSges(3,2) + Ifges(3,3) * t39 / 0.2e1) * t39 + (t9 * mrSges(4,1) - t10 * mrSges(4,2) + Ifges(4,3) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t40 * mrSges(2,3) + Ifges(2,5) * t49 + Ifges(2,1) * t43 / 0.2e1) * t43 + (t24 * mrSges(3,2) - t20 * mrSges(3,3) + Ifges(3,5) * t39 + Ifges(3,1) * t34 / 0.2e1) * t34 + (t18 * mrSges(4,2) - t9 * mrSges(4,3) + Ifges(4,5) * t31 + Ifges(4,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t41 * mrSges(2,3) + Ifges(2,4) * t43 + Ifges(2,6) * t49 + Ifges(2,2) * t42 / 0.2e1) * t42 + (-t24 * mrSges(3,1) + t21 * mrSges(3,3) + Ifges(3,4) * t34 + Ifges(3,6) * t39 + Ifges(3,2) * t33 / 0.2e1) * t33 + (-t18 * mrSges(4,1) + t10 * mrSges(4,3) + Ifges(4,4) * t27 + Ifges(4,6) * t31 + Ifges(4,2) * t26 / 0.2e1) * t26 + (t3 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t25) * t25 + (t7 * mrSges(5,2) + t5 * mrSges(6,2) - t3 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t23 + (Ifges(5,5) + Ifges(6,5)) * t25) * t23 + (-t7 * mrSges(5,1) - t5 * mrSges(6,1) + t4 * mrSges(5,3) + t2 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t22 + (Ifges(5,6) + Ifges(6,6)) * t25 + (Ifges(5,4) + Ifges(6,4)) * t23) * t22;
T = t6;

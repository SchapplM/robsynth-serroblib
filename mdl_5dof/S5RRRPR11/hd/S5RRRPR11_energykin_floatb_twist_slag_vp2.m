% Calculate kinetic energy for
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR11_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:20
% EndTime: 2019-12-31 21:32:21
% DurationCPUTime: 0.93s
% Computational Cost: add. (1339->129), mult. (1690->178), div. (0->0), fcn. (1236->8), ass. (0->47)
t60 = -pkin(3) - pkin(4);
t59 = cos(qJ(3));
t51 = sin(qJ(1));
t54 = cos(qJ(1));
t39 = -t51 * V_base(4) + t54 * V_base(5);
t40 = t51 * V_base(5) + t54 * V_base(4);
t25 = -pkin(1) * t39 - pkin(6) * t40 + V_base(3);
t45 = V_base(5) * pkin(5) + V_base(1);
t46 = -V_base(4) * pkin(5) + V_base(2);
t36 = t54 * t45 + t51 * t46;
t47 = V_base(6) + qJD(1);
t31 = pkin(6) * t47 + t36;
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t20 = t50 * t25 + t53 * t31;
t38 = qJD(2) - t39;
t17 = pkin(7) * t38 + t20;
t35 = -t51 * t45 + t46 * t54;
t30 = -pkin(1) * t47 - t35;
t33 = -t40 * t50 + t53 * t47;
t34 = t40 * t53 + t47 * t50;
t18 = -pkin(2) * t33 - pkin(7) * t34 + t30;
t49 = sin(qJ(3));
t10 = t59 * t17 + t49 * t18;
t19 = t53 * t25 - t50 * t31;
t32 = qJD(3) - t33;
t7 = t32 * qJ(4) + t10;
t58 = pkin(2) * t38 + t19;
t9 = -t49 * t17 + t59 * t18;
t57 = qJD(4) - t9;
t22 = t59 * t34 + t49 * t38;
t56 = qJ(4) * t22 + t58;
t55 = V_base(3) ^ 2;
t52 = cos(qJ(5));
t48 = sin(qJ(5));
t29 = qJD(5) - t32;
t21 = t34 * t49 - t59 * t38;
t12 = t21 * t48 + t22 * t52;
t11 = t21 * t52 - t22 * t48;
t8 = pkin(3) * t21 - t56;
t6 = -t32 * pkin(3) + t57;
t5 = t60 * t21 + t56;
t4 = pkin(8) * t21 + t7;
t3 = -t22 * pkin(8) + t60 * t32 + t57;
t2 = t3 * t48 + t4 * t52;
t1 = t3 * t52 - t4 * t48;
t13 = m(2) * (t35 ^ 2 + t36 ^ 2 + t55) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t55) / 0.2e1 + m(3) * (t19 ^ 2 + t20 ^ 2 + t30 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t58 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(5) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t35 * mrSges(2,1) - t36 * mrSges(2,2) + Ifges(2,3) * t47 / 0.2e1) * t47 + (t19 * mrSges(3,1) - t20 * mrSges(3,2) + Ifges(3,3) * t38 / 0.2e1) * t38 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t35 * mrSges(2,3) + Ifges(2,5) * t47 + Ifges(2,1) * t40 / 0.2e1) * t40 + (t30 * mrSges(3,2) - t19 * mrSges(3,3) + Ifges(3,5) * t38 + Ifges(3,1) * t34 / 0.2e1) * t34 + (t5 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t12 / 0.2e1) * t12 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t36 * mrSges(2,3) + Ifges(2,4) * t40 + Ifges(2,6) * t47 + Ifges(2,2) * t39 / 0.2e1) * t39 + (-t30 * mrSges(3,1) + t20 * mrSges(3,3) + Ifges(3,4) * t34 + Ifges(3,6) * t38 + Ifges(3,2) * t33 / 0.2e1) * t33 + (-t5 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t12 + Ifges(6,6) * t29 + Ifges(6,2) * t11 / 0.2e1) * t11 + (t9 * mrSges(4,1) - t6 * mrSges(5,1) - t10 * mrSges(4,2) + t7 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t32) * t32 + (-t58 * mrSges(4,2) + t6 * mrSges(5,2) - t9 * mrSges(4,3) - t8 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t22 + (Ifges(5,4) + Ifges(4,5)) * t32) * t22 + (-t58 * mrSges(4,1) + t8 * mrSges(5,1) - t7 * mrSges(5,2) - t10 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t21 + (-Ifges(4,6) + Ifges(5,6)) * t32 + (-Ifges(4,4) + Ifges(5,5)) * t22) * t21;
T = t13;

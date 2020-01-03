% Calculate kinetic energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR10_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR10_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:32
% EndTime: 2019-12-31 21:26:33
% DurationCPUTime: 1.11s
% Computational Cost: add. (3231->137), mult. (4866->203), div. (0->0), fcn. (4038->12), ass. (0->56)
t50 = V_base(5) * pkin(6) + V_base(1);
t51 = -V_base(4) * pkin(6) + V_base(2);
t60 = sin(qJ(1));
t64 = cos(qJ(1));
t43 = -t50 * t60 + t64 * t51;
t52 = V_base(6) + qJD(1);
t56 = cos(pkin(5));
t46 = t60 * V_base(5) + t64 * V_base(4);
t71 = pkin(7) * t46;
t38 = pkin(1) * t52 - t56 * t71 + t43;
t45 = -t60 * V_base(4) + t64 * V_base(5);
t54 = sin(pkin(5));
t41 = -pkin(1) * t45 - t54 * t71 + V_base(3);
t72 = t38 * t56 + t41 * t54;
t44 = t64 * t50 + t60 * t51;
t66 = t45 * t56 + t52 * t54;
t35 = t66 * pkin(7) + t44;
t59 = sin(qJ(2));
t63 = cos(qJ(2));
t26 = -t59 * t35 + t72 * t63;
t36 = -t46 * t59 + t66 * t63;
t28 = -t38 * t54 + t56 * t41;
t37 = t46 * t63 + t66 * t59;
t19 = -pkin(2) * t36 - pkin(8) * t37 + t28;
t27 = t63 * t35 + t72 * t59;
t42 = -t45 * t54 + t52 * t56 + qJD(2);
t22 = pkin(8) * t42 + t27;
t58 = sin(qJ(3));
t62 = cos(qJ(3));
t13 = t58 * t19 + t62 * t22;
t29 = -t37 * t58 + t42 * t62;
t11 = qJ(4) * t29 + t13;
t53 = sin(pkin(10));
t55 = cos(pkin(10));
t12 = t62 * t19 - t22 * t58;
t30 = t37 * t62 + t42 * t58;
t34 = qJD(3) - t36;
t9 = pkin(3) * t34 - qJ(4) * t30 + t12;
t6 = t55 * t11 + t53 * t9;
t5 = -t11 * t53 + t55 * t9;
t24 = t29 * t55 - t30 * t53;
t21 = -pkin(2) * t42 - t26;
t14 = -pkin(3) * t29 + qJD(4) + t21;
t65 = V_base(3) ^ 2;
t61 = cos(qJ(5));
t57 = sin(qJ(5));
t25 = t29 * t53 + t30 * t55;
t23 = qJD(5) - t24;
t16 = t25 * t61 + t34 * t57;
t15 = -t25 * t57 + t34 * t61;
t7 = -pkin(4) * t24 - pkin(9) * t25 + t14;
t4 = pkin(9) * t34 + t6;
t3 = -pkin(4) * t34 - t5;
t2 = t4 * t61 + t57 * t7;
t1 = -t4 * t57 + t61 * t7;
t8 = m(2) * (t43 ^ 2 + t44 ^ 2 + t65) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t65) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t13 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + Ifges(2,3) * t52 / 0.2e1) * t52 + (t26 * mrSges(3,1) - t27 * mrSges(3,2) + Ifges(3,3) * t42 / 0.2e1) * t42 + (t21 * mrSges(4,2) - t12 * mrSges(4,3) + Ifges(4,1) * t30 / 0.2e1) * t30 + (t14 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t25 / 0.2e1) * t25 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t23 / 0.2e1) * t23 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,5) * t52 + Ifges(2,1) * t46 / 0.2e1) * t46 + (t28 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,5) * t42 + Ifges(3,1) * t37 / 0.2e1) * t37 + (-t21 * mrSges(4,1) + t13 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,2) * t29 / 0.2e1) * t29 + (-t14 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t25 + Ifges(5,2) * t24 / 0.2e1) * t24 + (t3 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t23 + Ifges(6,1) * t16 / 0.2e1) * t16 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t46 + Ifges(2,6) * t52 + Ifges(2,2) * t45 / 0.2e1) * t45 + (-t28 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t37 + Ifges(3,6) * t42 + Ifges(3,2) * t36 / 0.2e1) * t36 + (-t3 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t16 + Ifges(6,6) * t23 + Ifges(6,2) * t15 / 0.2e1) * t15 + (t12 * mrSges(4,1) + t5 * mrSges(5,1) - t13 * mrSges(4,2) - t6 * mrSges(5,2) + Ifges(4,5) * t30 + Ifges(5,5) * t25 + Ifges(4,6) * t29 + Ifges(5,6) * t24 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t34) * t34;
T = t8;

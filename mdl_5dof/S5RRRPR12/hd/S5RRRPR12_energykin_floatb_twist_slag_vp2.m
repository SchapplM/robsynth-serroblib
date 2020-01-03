% Calculate kinetic energy for
% S5RRRPR12
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:36:58
% EndTime: 2019-12-31 21:36:59
% DurationCPUTime: 1.12s
% Computational Cost: add. (3201->137), mult. (4814->203), div. (0->0), fcn. (3992->12), ass. (0->56)
t52 = V_base(5) * pkin(6) + V_base(1);
t53 = -V_base(4) * pkin(6) + V_base(2);
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t45 = -t52 * t62 + t65 * t53;
t54 = V_base(6) + qJD(1);
t58 = cos(pkin(5));
t48 = t62 * V_base(5) + t65 * V_base(4);
t73 = pkin(7) * t48;
t39 = pkin(1) * t54 - t58 * t73 + t45;
t47 = -t62 * V_base(4) + t65 * V_base(5);
t56 = sin(pkin(5));
t42 = -pkin(1) * t47 - t56 * t73 + V_base(3);
t74 = t39 * t58 + t42 * t56;
t46 = t65 * t52 + t62 * t53;
t67 = t47 * t58 + t54 * t56;
t36 = t67 * pkin(7) + t46;
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t24 = -t61 * t36 + t74 * t64;
t37 = -t48 * t61 + t67 * t64;
t28 = -t39 * t56 + t58 * t42;
t38 = t48 * t64 + t67 * t61;
t19 = -pkin(2) * t37 - pkin(8) * t38 + t28;
t25 = t64 * t36 + t74 * t61;
t44 = -t47 * t56 + t54 * t58 + qJD(2);
t23 = pkin(8) * t44 + t25;
t60 = sin(qJ(3));
t72 = cos(qJ(3));
t12 = t60 * t19 + t72 * t23;
t35 = qJD(3) - t37;
t10 = qJ(4) * t35 + t12;
t22 = -pkin(2) * t44 - t24;
t30 = t38 * t60 - t72 * t44;
t31 = t72 * t38 + t60 * t44;
t15 = pkin(3) * t30 - qJ(4) * t31 + t22;
t55 = sin(pkin(10));
t57 = cos(pkin(10));
t6 = t57 * t10 + t55 * t15;
t5 = -t10 * t55 + t57 * t15;
t11 = t72 * t19 - t60 * t23;
t9 = -t35 * pkin(3) + qJD(4) - t11;
t66 = V_base(3) ^ 2;
t63 = cos(qJ(5));
t59 = sin(qJ(5));
t29 = qJD(5) + t30;
t27 = t31 * t57 + t35 * t55;
t26 = -t31 * t55 + t35 * t57;
t17 = t26 * t59 + t27 * t63;
t16 = t26 * t63 - t27 * t59;
t7 = -t26 * pkin(4) + t9;
t4 = pkin(9) * t26 + t6;
t3 = pkin(4) * t30 - pkin(9) * t27 + t5;
t2 = t3 * t59 + t4 * t63;
t1 = t3 * t63 - t4 * t59;
t8 = m(2) * (t45 ^ 2 + t46 ^ 2 + t66) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t66) / 0.2e1 + m(3) * (t24 ^ 2 + t25 ^ 2 + t28 ^ 2) / 0.2e1 + m(4) * (t11 ^ 2 + t12 ^ 2 + t22 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(5) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t45 * mrSges(2,1) - t46 * mrSges(2,2) + Ifges(2,3) * t54 / 0.2e1) * t54 + (t24 * mrSges(3,1) - t25 * mrSges(3,2) + Ifges(3,3) * t44 / 0.2e1) * t44 + (t11 * mrSges(4,1) - t12 * mrSges(4,2) + Ifges(4,3) * t35 / 0.2e1) * t35 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t29 / 0.2e1) * t29 + (t9 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,1) * t27 / 0.2e1) * t27 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t45 * mrSges(2,3) + Ifges(2,5) * t54 + Ifges(2,1) * t48 / 0.2e1) * t48 + (t28 * mrSges(3,2) - t24 * mrSges(3,3) + Ifges(3,5) * t44 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t22 * mrSges(4,2) - t11 * mrSges(4,3) + Ifges(4,5) * t35 + Ifges(4,1) * t31 / 0.2e1) * t31 + (-t9 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t27 + Ifges(5,2) * t26 / 0.2e1) * t26 + (t7 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t29 + Ifges(6,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t46 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t54 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t28 * mrSges(3,1) + t25 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t44 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t7 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t17 + Ifges(6,6) * t29 + Ifges(6,2) * t16 / 0.2e1) * t16 + (t22 * mrSges(4,1) + t5 * mrSges(5,1) - t6 * mrSges(5,2) - t12 * mrSges(4,3) - Ifges(4,4) * t31 + Ifges(5,5) * t27 - Ifges(4,6) * t35 + Ifges(5,6) * t26 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t30) * t30;
T = t8;

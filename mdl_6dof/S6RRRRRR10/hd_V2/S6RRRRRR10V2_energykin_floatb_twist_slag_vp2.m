% Calculate kinetic energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10V2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:26
% EndTime: 2019-04-11 14:41:27
% DurationCPUTime: 1.41s
% Computational Cost: add. (2745->146), mult. (3585->218), div. (0->0), fcn. (2936->12), ass. (0->57)
t58 = sin(qJ(1));
t64 = cos(qJ(1));
t47 = t58 * V_base(5) + t64 * V_base(4);
t52 = V_base(6) + qJD(1);
t57 = sin(qJ(2));
t63 = cos(qJ(2));
t38 = -t47 * t57 + t52 * t63;
t39 = t47 * t63 + t52 * t57;
t56 = sin(qJ(3));
t62 = cos(qJ(3));
t31 = t38 * t62 - t39 * t56;
t32 = t38 * t56 + t39 * t62;
t49 = V_base(5) * pkin(4) + V_base(1);
t50 = -V_base(4) * pkin(4) + V_base(2);
t40 = -t58 * t49 + t50 * t64;
t37 = -pkin(1) * t52 - t40;
t33 = -pkin(2) * t38 + t37;
t15 = -pkin(3) * t31 - pkin(5) * t32 + t33;
t41 = t49 * t64 + t50 * t58;
t46 = -t58 * V_base(4) + t64 * V_base(5);
t43 = -pkin(1) * t46 + V_base(3);
t35 = -t41 * t57 + t63 * t43;
t45 = qJD(2) - t46;
t29 = pkin(2) * t45 + t35;
t36 = t41 * t63 + t43 * t57;
t24 = t56 * t29 + t62 * t36;
t44 = qJD(3) + t45;
t22 = pkin(5) * t44 + t24;
t55 = sin(qJ(4));
t61 = cos(qJ(4));
t11 = t15 * t55 + t22 * t61;
t23 = t29 * t62 - t36 * t56;
t21 = -pkin(3) * t44 - t23;
t54 = sin(qJ(5));
t60 = cos(qJ(5));
t5 = t11 * t54 - t60 * t21;
t67 = t5 ^ 2;
t9 = -t61 * t15 + t22 * t55;
t66 = t9 ^ 2;
t7 = t60 * t11 + t54 * t21;
t27 = t32 * t61 + t44 * t55;
t30 = qJD(4) - t31;
t17 = -t27 * t54 + t30 * t60;
t26 = -t32 * t55 + t44 * t61;
t65 = V_base(3) ^ 2;
t59 = cos(qJ(6));
t53 = sin(qJ(6));
t25 = qJD(5) - t26;
t18 = t27 * t60 + t30 * t54;
t16 = qJD(6) - t17;
t13 = t18 * t59 + t25 * t53;
t12 = -t18 * t53 + t25 * t59;
t4 = -pkin(6) * t18 + t9;
t3 = pkin(6) * t25 + t7;
t2 = t3 * t59 + t4 * t53;
t1 = -t3 * t53 + t4 * t59;
t6 = (-t9 * mrSges(6,1) + t7 * mrSges(6,3) + Ifges(6,4) * t18 + Ifges(6,6) * t25 + Ifges(6,2) * t17 / 0.2e1) * t17 + m(5) * (t11 ^ 2 + t21 ^ 2 + t66) / 0.2e1 + m(6) * (t7 ^ 2 + t66 + t67) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t67) / 0.2e1 + (-t5 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t13 + Ifges(7,6) * t16 + Ifges(7,2) * t12 / 0.2e1) * t12 + (t5 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t16 + Ifges(7,1) * t13 / 0.2e1) * t13 + (t23 * mrSges(4,1) - t24 * mrSges(4,2) + Ifges(4,3) * t44 / 0.2e1) * t44 + (t40 * mrSges(2,1) - t41 * mrSges(2,2) + Ifges(2,3) * t52 / 0.2e1) * t52 + (t35 * mrSges(3,1) - t36 * mrSges(3,2) + Ifges(3,3) * t45 / 0.2e1) * t45 + (t33 * mrSges(4,2) - t23 * mrSges(4,3) + Ifges(4,5) * t44 + Ifges(4,1) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(2,1) + t41 * mrSges(2,3) + Ifges(2,4) * t47 + Ifges(2,6) * t52 + Ifges(2,2) * t46 / 0.2e1) * t46 + (-t33 * mrSges(4,1) + t24 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t44 + Ifges(4,2) * t31 / 0.2e1) * t31 + (V_base(3) * mrSges(2,2) - t40 * mrSges(2,3) + Ifges(2,5) * t52 + Ifges(2,1) * t47 / 0.2e1) * t47 + (t37 * mrSges(3,2) - t35 * mrSges(3,3) + Ifges(3,5) * t45 + Ifges(3,1) * t39 / 0.2e1) * t39 + (-t9 * mrSges(5,1) - t11 * mrSges(5,2) + Ifges(5,3) * t30 / 0.2e1) * t30 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (-t5 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,3) * t25 / 0.2e1) * t25 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + m(2) * (t40 ^ 2 + t41 ^ 2 + t65) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t65) / 0.2e1 + m(4) * (t23 ^ 2 + t24 ^ 2 + t33 ^ 2) / 0.2e1 + m(3) * (t35 ^ 2 + t36 ^ 2 + t37 ^ 2) / 0.2e1 + (-t21 * mrSges(5,1) + t11 * mrSges(5,3) + Ifges(5,4) * t27 + Ifges(5,6) * t30 + Ifges(5,2) * t26 / 0.2e1) * t26 + (t21 * mrSges(5,2) + t9 * mrSges(5,3) + Ifges(5,5) * t30 + Ifges(5,1) * t27 / 0.2e1) * t27 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t16 / 0.2e1) * t16 + (-t37 * mrSges(3,1) + t36 * mrSges(3,3) + Ifges(3,4) * t39 + Ifges(3,6) * t45 + Ifges(3,2) * t38 / 0.2e1) * t38 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t9 * mrSges(6,2) + t5 * mrSges(6,3) + t25 * Ifges(6,5) + Ifges(6,1) * t18 / 0.2e1) * t18;
T  = t6;

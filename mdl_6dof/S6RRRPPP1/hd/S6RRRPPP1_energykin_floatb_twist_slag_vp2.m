% Calculate kinetic energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPP1_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:23
% EndTime: 2019-03-09 15:14:24
% DurationCPUTime: 0.99s
% Computational Cost: add. (3413->156), mult. (4469->210), div. (0->0), fcn. (3618->10), ass. (0->57)
t72 = pkin(4) + qJ(6);
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t47 = -t60 * V_base(4) + t63 * V_base(5);
t48 = t60 * V_base(5) + t63 * V_base(4);
t36 = -pkin(1) * t47 - pkin(8) * t48 + V_base(3);
t52 = V_base(5) * pkin(7) + V_base(1);
t53 = -V_base(4) * pkin(7) + V_base(2);
t45 = t63 * t52 + t60 * t53;
t54 = V_base(6) + qJD(1);
t40 = pkin(8) * t54 + t45;
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t30 = t59 * t36 + t62 * t40;
t46 = qJD(2) - t47;
t26 = pkin(9) * t46 + t30;
t44 = -t60 * t52 + t53 * t63;
t39 = -pkin(1) * t54 - t44;
t42 = -t48 * t59 + t54 * t62;
t43 = t48 * t62 + t54 * t59;
t28 = -pkin(2) * t42 - pkin(9) * t43 + t39;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t19 = t61 * t26 + t58 * t28;
t33 = t43 * t61 + t46 * t58;
t71 = qJ(4) * t33;
t70 = cos(pkin(10));
t32 = -t43 * t58 + t46 * t61;
t41 = qJD(3) - t42;
t56 = sin(pkin(6));
t57 = cos(pkin(6));
t67 = t32 * t57 + t41 * t56;
t13 = qJ(4) * t67 + t19;
t18 = -t26 * t58 + t61 * t28;
t14 = pkin(3) * t41 - t57 * t71 + t18;
t29 = t36 * t62 - t59 * t40;
t25 = -pkin(2) * t46 - t29;
t17 = -pkin(3) * t32 - t56 * t71 + t25;
t55 = sin(pkin(10));
t8 = t70 * t13 + (t14 * t57 + t17 * t56) * t55;
t69 = t56 * t70;
t68 = t57 * t70;
t9 = -t14 * t56 + t57 * t17 + qJD(4);
t27 = -t32 * t56 + t41 * t57;
t5 = -qJ(5) * t27 - t8;
t21 = t33 * t70 + t55 * t67;
t66 = -qJ(5) * t21 + t9;
t7 = -t55 * t13 + t14 * t68 + t17 * t69;
t65 = qJD(5) - t7;
t64 = V_base(3) ^ 2;
t20 = -t32 * t68 + t33 * t55 - t41 * t69;
t6 = pkin(4) * t20 + t66;
t4 = -t27 * pkin(4) + t65;
t3 = t20 * t72 + t66;
t2 = -pkin(5) * t20 + qJD(6) - t5;
t1 = t21 * pkin(5) - t72 * t27 + t65;
t10 = m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t64) / 0.2e1 + m(2) * (t44 ^ 2 + t45 ^ 2 + t64) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t39 ^ 2) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t25 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t7 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + Ifges(2,3) * t54 / 0.2e1) * t54 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t46 / 0.2e1) * t46 + (t18 * mrSges(4,1) - t19 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,5) * t54 + Ifges(2,1) * t48 / 0.2e1) * t48 + (t39 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t46 + Ifges(3,1) * t43 / 0.2e1) * t43 + (t25 * mrSges(4,2) - t18 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t33 / 0.2e1) * t33 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t54 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t39 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t43 + Ifges(3,6) * t46 + Ifges(3,2) * t42 / 0.2e1) * t42 + (-t25 * mrSges(4,1) + t19 * mrSges(4,3) + Ifges(4,4) * t33 + Ifges(4,6) * t41 + Ifges(4,2) * t32 / 0.2e1) * t32 + (t7 * mrSges(5,1) - t8 * mrSges(5,2) + t4 * mrSges(6,2) + t2 * mrSges(7,2) - t5 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(5,3) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t27) * t27 + (t4 * mrSges(6,1) + t1 * mrSges(7,1) + t9 * mrSges(5,2) - t3 * mrSges(7,2) - t7 * mrSges(5,3) - t6 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t21 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t27) * t21 + (t9 * mrSges(5,1) + t5 * mrSges(6,1) - t2 * mrSges(7,1) - t6 * mrSges(6,2) - t8 * mrSges(5,3) + t3 * mrSges(7,3) + (Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t20 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t27 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t21) * t20;
T  = t10;

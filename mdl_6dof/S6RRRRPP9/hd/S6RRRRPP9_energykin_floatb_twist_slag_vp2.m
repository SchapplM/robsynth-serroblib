% Calculate kinetic energy for
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:51
% EndTime: 2019-03-09 21:41:52
% DurationCPUTime: 0.99s
% Computational Cost: add. (3539->154), mult. (5277->206), div. (0->0), fcn. (4376->10), ass. (0->56)
t51 = V_base(5) * pkin(7) + V_base(1);
t52 = -V_base(4) * pkin(7) + V_base(2);
t59 = sin(qJ(1));
t62 = cos(qJ(1));
t44 = -t51 * t59 + t62 * t52;
t53 = V_base(6) + qJD(1);
t55 = cos(pkin(6));
t47 = t59 * V_base(5) + t62 * V_base(4);
t73 = pkin(8) * t47;
t39 = pkin(1) * t53 - t55 * t73 + t44;
t46 = -t59 * V_base(4) + t62 * V_base(5);
t54 = sin(pkin(6));
t42 = -pkin(1) * t46 - t54 * t73 + V_base(3);
t74 = t39 * t55 + t42 * t54;
t45 = t62 * t51 + t59 * t52;
t66 = t46 * t55 + t53 * t54;
t36 = pkin(8) * t66 + t45;
t58 = sin(qJ(2));
t61 = cos(qJ(2));
t23 = -t58 * t36 + t61 * t74;
t37 = -t47 * t58 + t61 * t66;
t72 = cos(qJ(4));
t68 = pkin(4) + qJ(6);
t27 = -t39 * t54 + t55 * t42;
t38 = t47 * t61 + t58 * t66;
t18 = -pkin(2) * t37 - pkin(9) * t38 + t27;
t24 = t61 * t36 + t74 * t58;
t43 = -t46 * t54 + t53 * t55 + qJD(2);
t22 = pkin(9) * t43 + t24;
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t14 = t57 * t18 + t60 * t22;
t35 = qJD(3) - t37;
t12 = pkin(10) * t35 + t14;
t21 = -pkin(2) * t43 - t23;
t29 = -t38 * t57 + t43 * t60;
t30 = t38 * t60 + t43 * t57;
t16 = -pkin(3) * t29 - pkin(10) * t30 + t21;
t56 = sin(qJ(4));
t7 = t72 * t12 + t56 * t16;
t13 = t18 * t60 - t57 * t22;
t28 = qJD(4) - t29;
t4 = -qJ(5) * t28 - t7;
t6 = -t56 * t12 + t16 * t72;
t11 = -pkin(3) * t35 - t13;
t65 = qJD(5) - t6;
t26 = t30 * t72 + t56 * t35;
t64 = -qJ(5) * t26 + t11;
t63 = V_base(3) ^ 2;
t25 = t30 * t56 - t35 * t72;
t8 = pkin(4) * t25 + t64;
t5 = t25 * t68 + t64;
t3 = -t28 * pkin(4) + t65;
t2 = -pkin(5) * t25 + qJD(6) - t4;
t1 = t26 * pkin(5) - t28 * t68 + t65;
t9 = m(2) * (t44 ^ 2 + t45 ^ 2 + t63) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t63) / 0.2e1 + m(3) * (t23 ^ 2 + t24 ^ 2 + t27 ^ 2) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t21 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + Ifges(2,3) * t53 / 0.2e1) * t53 + (t23 * mrSges(3,1) - t24 * mrSges(3,2) + Ifges(3,3) * t43 / 0.2e1) * t43 + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t35 / 0.2e1) * t35 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,5) * t53 + Ifges(2,1) * t47 / 0.2e1) * t47 + (t27 * mrSges(3,2) - t23 * mrSges(3,3) + Ifges(3,5) * t43 + Ifges(3,1) * t38 / 0.2e1) * t38 + (t21 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t35 + Ifges(4,1) * t30 / 0.2e1) * t30 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t47 + Ifges(2,6) * t53 + Ifges(2,2) * t46 / 0.2e1) * t46 + (-t27 * mrSges(3,1) + t24 * mrSges(3,3) + Ifges(3,4) * t38 + Ifges(3,6) * t43 + Ifges(3,2) * t37 / 0.2e1) * t37 + (-t21 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t30 + Ifges(4,6) * t35 + Ifges(4,2) * t29 / 0.2e1) * t29 + (t6 * mrSges(5,1) - t7 * mrSges(5,2) + t3 * mrSges(6,2) + t2 * mrSges(7,2) - t4 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t28) * t28 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) + t11 * mrSges(5,2) - t5 * mrSges(7,2) - t6 * mrSges(5,3) - t8 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t26 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t28) * t26 + (t11 * mrSges(5,1) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t8 * mrSges(6,2) - t7 * mrSges(5,3) + t5 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t25 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t28 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t26) * t25;
T  = t9;

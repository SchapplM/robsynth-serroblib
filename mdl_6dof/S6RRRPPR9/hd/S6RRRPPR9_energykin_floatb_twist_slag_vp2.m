% Calculate kinetic energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:31
% EndTime: 2019-03-09 16:11:32
% DurationCPUTime: 1.36s
% Computational Cost: add. (4301->158), mult. (6449->221), div. (0->0), fcn. (5406->12), ass. (0->62)
t58 = V_base(5) * pkin(7) + V_base(1);
t59 = -V_base(4) * pkin(7) + V_base(2);
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t51 = -t58 * t67 + t70 * t59;
t60 = V_base(6) + qJD(1);
t63 = cos(pkin(6));
t54 = t67 * V_base(5) + t70 * V_base(4);
t82 = pkin(8) * t54;
t45 = pkin(1) * t60 - t63 * t82 + t51;
t53 = -t67 * V_base(4) + t70 * V_base(5);
t62 = sin(pkin(6));
t48 = -pkin(1) * t53 - t62 * t82 + V_base(3);
t84 = t45 * t63 + t48 * t62;
t52 = t70 * t58 + t67 * t59;
t75 = t53 * t63 + t60 * t62;
t42 = pkin(8) * t75 + t52;
t66 = sin(qJ(2));
t69 = cos(qJ(2));
t28 = -t66 * t42 + t69 * t84;
t43 = -t54 * t66 + t69 * t75;
t83 = -pkin(4) - pkin(5);
t81 = cos(qJ(3));
t32 = -t45 * t62 + t63 * t48;
t44 = t54 * t69 + t66 * t75;
t23 = -pkin(2) * t43 - pkin(9) * t44 + t32;
t29 = t69 * t42 + t66 * t84;
t50 = -t53 * t62 + t60 * t63 + qJD(2);
t27 = pkin(9) * t50 + t29;
t65 = sin(qJ(3));
t16 = t65 * t23 + t81 * t27;
t41 = qJD(3) - t43;
t14 = qJ(4) * t41 + t16;
t26 = -pkin(2) * t50 - t28;
t35 = t44 * t65 - t81 * t50;
t36 = t44 * t81 + t65 * t50;
t18 = pkin(3) * t35 - qJ(4) * t36 + t26;
t61 = sin(pkin(11));
t77 = cos(pkin(11));
t9 = t77 * t14 + t61 * t18;
t15 = t81 * t23 - t65 * t27;
t6 = t35 * qJ(5) + t9;
t8 = -t61 * t14 + t18 * t77;
t74 = pkin(3) * t41 - qJD(4) + t15;
t73 = qJD(5) - t8;
t31 = t36 * t77 + t61 * t41;
t72 = qJ(5) * t31 + t74;
t71 = V_base(3) ^ 2;
t68 = cos(qJ(6));
t64 = sin(qJ(6));
t34 = qJD(6) - t35;
t30 = t36 * t61 - t41 * t77;
t20 = t30 * t64 + t31 * t68;
t19 = t30 * t68 - t31 * t64;
t10 = pkin(4) * t30 - t72;
t7 = t30 * t83 + t72;
t5 = -t35 * pkin(4) + t73;
t4 = pkin(10) * t30 + t6;
t3 = -t31 * pkin(10) + t35 * t83 + t73;
t2 = t3 * t64 + t4 * t68;
t1 = t3 * t68 - t4 * t64;
t11 = m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t71) / 0.2e1 + m(2) * (t51 ^ 2 + t52 ^ 2 + t71) / 0.2e1 + m(3) * (t28 ^ 2 + t29 ^ 2 + t32 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t10 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t51 * mrSges(2,1) - t52 * mrSges(2,2) + Ifges(2,3) * t60 / 0.2e1) * t60 + (t28 * mrSges(3,1) - t29 * mrSges(3,2) + Ifges(3,3) * t50 / 0.2e1) * t50 + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t41 / 0.2e1) * t41 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t51 * mrSges(2,3) + Ifges(2,5) * t60 + Ifges(2,1) * t54 / 0.2e1) * t54 + (t32 * mrSges(3,2) - t28 * mrSges(3,3) + Ifges(3,5) * t50 + Ifges(3,1) * t44 / 0.2e1) * t44 + (t26 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,5) * t41 + Ifges(4,1) * t36 / 0.2e1) * t36 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t34 + Ifges(7,1) * t20 / 0.2e1) * t20 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t52 * mrSges(2,3) + Ifges(2,4) * t54 + Ifges(2,6) * t60 + Ifges(2,2) * t53 / 0.2e1) * t53 + (-t32 * mrSges(3,1) + t29 * mrSges(3,3) + Ifges(3,4) * t44 + Ifges(3,6) * t50 + Ifges(3,2) * t43 / 0.2e1) * t43 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t20 + Ifges(7,6) * t34 + Ifges(7,2) * t19 / 0.2e1) * t19 + (-t74 * mrSges(5,2) + t5 * mrSges(6,2) - t8 * mrSges(5,3) - t10 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t31) * t31 + (-t74 * mrSges(5,1) + t10 * mrSges(6,1) - t6 * mrSges(6,2) - t9 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t30 + (-Ifges(5,4) + Ifges(6,5)) * t31) * t30 + (t26 * mrSges(4,1) + t8 * mrSges(5,1) - t5 * mrSges(6,1) - t9 * mrSges(5,2) - t16 * mrSges(4,3) + t6 * mrSges(6,3) - Ifges(4,4) * t36 - Ifges(4,6) * t41 + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t35 + (Ifges(6,4) + Ifges(5,5)) * t31 + (-Ifges(5,6) + Ifges(6,6)) * t30) * t35;
T  = t11;

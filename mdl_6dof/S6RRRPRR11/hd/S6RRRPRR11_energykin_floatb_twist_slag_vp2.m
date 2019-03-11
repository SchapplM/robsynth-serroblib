% Calculate kinetic energy for
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR11_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:43
% EndTime: 2019-03-09 19:24:44
% DurationCPUTime: 1.28s
% Computational Cost: add. (4029->160), mult. (6001->227), div. (0->0), fcn. (5006->12), ass. (0->62)
t72 = sin(qJ(1));
t76 = cos(qJ(1));
t58 = t72 * V_base(5) + t76 * V_base(4);
t83 = pkin(8) * t58;
t63 = pkin(7) * V_base(5) + V_base(1);
t64 = -pkin(7) * V_base(4) + V_base(2);
t54 = -t63 * t72 + t64 * t76;
t65 = V_base(6) + qJD(1);
t67 = cos(pkin(6));
t46 = pkin(1) * t65 - t67 * t83 + t54;
t57 = -t72 * V_base(4) + t76 * V_base(5);
t66 = sin(pkin(6));
t50 = -pkin(1) * t57 - t66 * t83 + V_base(3);
t31 = -t46 * t66 + t67 * t50;
t71 = sin(qJ(2));
t75 = cos(qJ(2));
t80 = t67 * t75;
t81 = t66 * t75;
t44 = t57 * t80 - t58 * t71 + t65 * t81;
t79 = t57 * t67 + t65 * t66;
t45 = t58 * t75 + t71 * t79;
t21 = -pkin(2) * t44 - pkin(9) * t45 + t31;
t55 = t63 * t76 + t64 * t72;
t43 = pkin(8) * t79 + t55;
t30 = t75 * t43 + (t46 * t67 + t50 * t66) * t71;
t53 = -t57 * t66 + t65 * t67 + qJD(2);
t25 = pkin(9) * t53 + t30;
t70 = sin(qJ(3));
t82 = cos(qJ(3));
t16 = t21 * t70 + t25 * t82;
t42 = qJD(3) - t44;
t14 = t42 * qJ(4) + t16;
t34 = t45 * t70 - t53 * t82;
t11 = pkin(10) * t34 + t14;
t69 = sin(qJ(5));
t74 = cos(qJ(5));
t35 = t45 * t82 + t53 * t70;
t15 = t21 * t82 - t25 * t70;
t78 = qJD(4) - t15;
t9 = -t35 * pkin(10) + (-pkin(3) - pkin(4)) * t42 + t78;
t6 = t11 * t74 + t69 * t9;
t29 = -t43 * t71 + t46 * t80 + t50 * t81;
t24 = -pkin(2) * t53 - t29;
t5 = -t11 * t69 + t74 * t9;
t27 = t34 * t74 - t35 * t69;
t17 = pkin(3) * t34 - qJ(4) * t35 + t24;
t12 = -pkin(4) * t34 - t17;
t77 = V_base(3) ^ 2;
t73 = cos(qJ(6));
t68 = sin(qJ(6));
t41 = qJD(5) - t42;
t28 = t34 * t69 + t35 * t74;
t26 = qJD(6) - t27;
t19 = t28 * t73 + t41 * t68;
t18 = -t28 * t68 + t41 * t73;
t13 = -pkin(3) * t42 + t78;
t7 = -pkin(5) * t27 - pkin(11) * t28 + t12;
t4 = pkin(11) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t4 * t73 + t68 * t7;
t1 = -t4 * t68 + t7 * t73;
t8 = (V_base(3) * mrSges(2,2) - t54 * mrSges(2,3) + Ifges(2,5) * t65 + Ifges(2,1) * t58 / 0.2e1) * t58 + m(2) * (t54 ^ 2 + t55 ^ 2 + t77) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t77) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t24 ^ 2) / 0.2e1 + m(5) * (t13 ^ 2 + t14 ^ 2 + t17 ^ 2) / 0.2e1 + m(6) * (t12 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t55 * mrSges(2,3) + Ifges(2,4) * t58 + Ifges(2,6) * t65 + Ifges(2,2) * t57 / 0.2e1) * t57 + (t31 * mrSges(3,2) - t29 * mrSges(3,3) + Ifges(3,5) * t53 + Ifges(3,1) * t45 / 0.2e1) * t45 + (t24 * mrSges(4,1) + t17 * mrSges(5,1) - t14 * mrSges(5,2) - t16 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t34 + (-Ifges(4,6) + Ifges(5,6)) * t42 + (-Ifges(4,4) + Ifges(5,5)) * t35) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t31 * mrSges(3,1) + t30 * mrSges(3,3) + Ifges(3,4) * t45 + Ifges(3,6) * t53 + Ifges(3,2) * t44 / 0.2e1) * t44 + (t24 * mrSges(4,2) + t13 * mrSges(5,2) - t15 * mrSges(4,3) - t17 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t35 + (Ifges(5,4) + Ifges(4,5)) * t42) * t35 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t41 / 0.2e1) * t41 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t26 + Ifges(7,1) * t19 / 0.2e1) * t19 + (-t12 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t41 + Ifges(6,2) * t27 / 0.2e1) * t27 + (t15 * mrSges(4,1) - t13 * mrSges(5,1) - t16 * mrSges(4,2) + t14 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t42) * t42 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t26 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t54 * mrSges(2,1) - t55 * mrSges(2,2) + Ifges(2,3) * t65 / 0.2e1) * t65 + (t12 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t41 + Ifges(6,1) * t28 / 0.2e1) * t28 + (t29 * mrSges(3,1) - t30 * mrSges(3,2) + Ifges(3,3) * t53 / 0.2e1) * t53 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t26 / 0.2e1) * t26;
T  = t8;

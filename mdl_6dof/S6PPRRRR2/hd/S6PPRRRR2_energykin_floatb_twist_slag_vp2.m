% Calculate kinetic energy for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:09
% EndTime: 2019-03-08 19:03:11
% DurationCPUTime: 1.79s
% Computational Cost: add. (8567->166), mult. (16221->249), div. (0->0), fcn. (14174->16), ass. (0->72)
t68 = qJ(1) * V_base(5) + V_base(1);
t69 = -qJ(1) * V_base(4) + V_base(2);
t72 = sin(pkin(12));
t76 = cos(pkin(12));
t61 = -t68 * t72 + t76 * t69;
t78 = cos(pkin(6));
t64 = t72 * V_base(5) + t76 * V_base(4);
t90 = qJ(2) * t64;
t55 = V_base(6) * pkin(1) - t78 * t90 + t61;
t63 = -t72 * V_base(4) + t76 * V_base(5);
t70 = V_base(3) + qJD(1);
t74 = sin(pkin(6));
t59 = -pkin(1) * t63 - t74 * t90 + t70;
t98 = t55 * t78 + t59 * t74;
t62 = t76 * t68 + t72 * t69;
t87 = t63 * t78 + t74 * V_base(6);
t52 = qJ(2) * t87 + t62;
t71 = sin(pkin(13));
t75 = cos(pkin(13));
t41 = -t52 * t71 + t98 * t75;
t60 = -t63 * t74 + t78 * V_base(6);
t77 = cos(pkin(7));
t54 = t64 * t75 + t71 * t87;
t96 = pkin(8) * t54;
t34 = pkin(2) * t60 - t77 * t96 + t41;
t47 = -t55 * t74 + t78 * t59 + qJD(2);
t53 = -t64 * t71 + t75 * t87;
t73 = sin(pkin(7));
t40 = -pkin(2) * t53 - t73 * t96 + t47;
t97 = t34 * t77 + t40 * t73;
t42 = t75 * t52 + t98 * t71;
t88 = t53 * t77 + t60 * t73;
t32 = pkin(8) * t88 + t42;
t82 = sin(qJ(3));
t86 = cos(qJ(3));
t24 = -t82 * t32 + t86 * t97;
t45 = -t82 * t54 + t86 * t88;
t25 = t86 * t32 + t97 * t82;
t48 = -t53 * t73 + t60 * t77 + qJD(3);
t19 = pkin(9) * t48 + t25;
t26 = -t34 * t73 + t77 * t40;
t46 = t54 * t86 + t82 * t88;
t23 = -pkin(3) * t45 - pkin(9) * t46 + t26;
t81 = sin(qJ(4));
t85 = cos(qJ(4));
t12 = t85 * t19 + t81 * t23;
t44 = qJD(4) - t45;
t10 = pkin(10) * t44 + t12;
t18 = -t48 * pkin(3) - t24;
t36 = -t81 * t46 + t48 * t85;
t37 = t46 * t85 + t48 * t81;
t15 = -t36 * pkin(4) - t37 * pkin(10) + t18;
t80 = sin(qJ(5));
t84 = cos(qJ(5));
t6 = t84 * t10 + t80 * t15;
t5 = -t10 * t80 + t84 * t15;
t11 = -t81 * t19 + t23 * t85;
t35 = qJD(5) - t36;
t9 = -pkin(4) * t44 - t11;
t83 = cos(qJ(6));
t79 = sin(qJ(6));
t33 = qJD(6) + t35;
t28 = t37 * t84 + t44 * t80;
t27 = -t37 * t80 + t44 * t84;
t21 = t27 * t79 + t28 * t83;
t20 = t27 * t83 - t28 * t79;
t7 = -pkin(5) * t27 + t9;
t4 = pkin(11) * t27 + t6;
t3 = pkin(5) * t35 - pkin(11) * t28 + t5;
t2 = t3 * t79 + t4 * t83;
t1 = t3 * t83 - t4 * t79;
t8 = m(3) * (t41 ^ 2 + t42 ^ 2 + t47 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t18 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (-t26 * mrSges(4,1) + t25 * mrSges(4,3) + Ifges(4,4) * t46 + Ifges(4,6) * t48 + Ifges(4,2) * t45 / 0.2e1) * t45 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t33 + Ifges(7,1) * t21 / 0.2e1) * t21 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t35 + Ifges(6,1) * t28 / 0.2e1) * t28 + (-t18 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t37 + Ifges(5,6) * t44 + Ifges(5,2) * t36 / 0.2e1) * t36 + (-t70 * mrSges(2,1) + t62 * mrSges(2,3) + Ifges(2,4) * t64 + Ifges(2,2) * t63 / 0.2e1) * t63 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t35 + Ifges(6,2) * t27 / 0.2e1) * t27 + (t47 * mrSges(3,2) - t41 * mrSges(3,3) + Ifges(3,5) * t60 + Ifges(3,1) * t54 / 0.2e1) * t54 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t21 + Ifges(7,6) * t33 + Ifges(7,2) * t20 / 0.2e1) * t20 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t35 / 0.2e1) * t35 + (t18 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t44 + Ifges(5,1) * t37 / 0.2e1) * t37 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t26 * mrSges(4,2) - t24 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t46 / 0.2e1) * t46 + (t24 * mrSges(4,1) - t25 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t33 / 0.2e1) * t33 + (-t47 * mrSges(3,1) + t42 * mrSges(3,3) + Ifges(3,4) * t54 + Ifges(3,6) * t60 + Ifges(3,2) * t53 / 0.2e1) * t53 + (t41 * mrSges(3,1) - t42 * mrSges(3,2) + Ifges(3,3) * t60 / 0.2e1) * t60 + (V_base(2) * mrSges(1,1) + t61 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t62 * mrSges(2,2) + Ifges(1,5) * V_base(4) + Ifges(2,5) * t64 + Ifges(1,6) * V_base(5) + Ifges(2,6) * t63 + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t61 ^ 2 + t62 ^ 2 + t70 ^ 2) / 0.2e1 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t44 / 0.2e1) * t44 + (t70 * mrSges(2,2) - t61 * mrSges(2,3) + Ifges(2,1) * t64 / 0.2e1) * t64;
T  = t8;

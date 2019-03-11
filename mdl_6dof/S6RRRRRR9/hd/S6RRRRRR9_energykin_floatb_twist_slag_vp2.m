% Calculate kinetic energy for
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR9_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:36
% EndTime: 2019-03-10 05:16:38
% DurationCPUTime: 1.91s
% Computational Cost: add. (10335->166), mult. (16221->250), div. (0->0), fcn. (14174->16), ass. (0->73)
t68 = pkin(8) * V_base(5) + V_base(1);
t69 = -pkin(8) * V_base(4) + V_base(2);
t80 = sin(qJ(1));
t86 = cos(qJ(1));
t61 = -t68 * t80 + t69 * t86;
t70 = V_base(6) + qJD(1);
t74 = cos(pkin(6));
t64 = t80 * V_base(5) + t86 * V_base(4);
t97 = pkin(9) * t64;
t55 = pkin(1) * t70 - t74 * t97 + t61;
t63 = -t80 * V_base(4) + t86 * V_base(5);
t72 = sin(pkin(6));
t59 = -pkin(1) * t63 - t72 * t97 + V_base(3);
t99 = t55 * t74 + t59 * t72;
t62 = t68 * t86 + t69 * t80;
t88 = t63 * t74 + t70 * t72;
t52 = pkin(9) * t88 + t62;
t79 = sin(qJ(2));
t85 = cos(qJ(2));
t41 = -t52 * t79 + t85 * t99;
t60 = -t63 * t72 + t70 * t74 + qJD(2);
t73 = cos(pkin(7));
t54 = t64 * t85 + t79 * t88;
t96 = pkin(10) * t54;
t37 = pkin(2) * t60 - t73 * t96 + t41;
t47 = -t55 * t72 + t59 * t74;
t53 = -t64 * t79 + t85 * t88;
t71 = sin(pkin(7));
t40 = -pkin(2) * t53 - t71 * t96 + t47;
t98 = t37 * t73 + t40 * t71;
t42 = t85 * t52 + t79 * t99;
t89 = t53 * t73 + t60 * t71;
t36 = pkin(10) * t89 + t42;
t78 = sin(qJ(3));
t84 = cos(qJ(3));
t24 = -t36 * t78 + t84 * t98;
t45 = -t54 * t78 + t84 * t89;
t25 = t84 * t36 + t78 * t98;
t48 = -t53 * t71 + t60 * t73 + qJD(3);
t21 = pkin(11) * t48 + t25;
t26 = -t37 * t71 + t40 * t73;
t46 = t54 * t84 + t78 * t89;
t23 = -pkin(3) * t45 - pkin(11) * t46 + t26;
t77 = sin(qJ(4));
t83 = cos(qJ(4));
t12 = t21 * t83 + t23 * t77;
t44 = qJD(4) - t45;
t10 = pkin(12) * t44 + t12;
t20 = -pkin(3) * t48 - t24;
t34 = -t46 * t77 + t48 * t83;
t35 = t46 * t83 + t48 * t77;
t15 = -pkin(4) * t34 - pkin(12) * t35 + t20;
t76 = sin(qJ(5));
t82 = cos(qJ(5));
t6 = t10 * t82 + t15 * t76;
t5 = -t10 * t76 + t15 * t82;
t11 = -t21 * t77 + t23 * t83;
t33 = qJD(5) - t34;
t9 = -pkin(4) * t44 - t11;
t87 = V_base(3) ^ 2;
t81 = cos(qJ(6));
t75 = sin(qJ(6));
t29 = qJD(6) + t33;
t28 = t35 * t82 + t44 * t76;
t27 = -t35 * t76 + t44 * t82;
t17 = t27 * t75 + t28 * t81;
t16 = t27 * t81 - t28 * t75;
t7 = -pkin(5) * t27 + t9;
t4 = pkin(13) * t27 + t6;
t3 = pkin(5) * t33 - pkin(13) * t28 + t5;
t2 = t3 * t75 + t4 * t81;
t1 = t3 * t81 - t4 * t75;
t8 = (-t47 * mrSges(3,1) + t42 * mrSges(3,3) + Ifges(3,4) * t54 + Ifges(3,6) * t60 + Ifges(3,2) * t53 / 0.2e1) * t53 + (t61 * mrSges(2,1) - t62 * mrSges(2,2) + Ifges(2,3) * t70 / 0.2e1) * t70 + (t24 * mrSges(4,1) - t25 * mrSges(4,2) + Ifges(4,3) * t48 / 0.2e1) * t48 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,6) * t33 + Ifges(6,2) * t27 / 0.2e1) * t27 + (t20 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t44 + Ifges(5,1) * t35 / 0.2e1) * t35 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t29 / 0.2e1) * t29 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t33 / 0.2e1) * t33 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t44 / 0.2e1) * t44 + (t41 * mrSges(3,1) - t42 * mrSges(3,2) + Ifges(3,3) * t60 / 0.2e1) * t60 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t29 + Ifges(7,2) * t16 / 0.2e1) * t16 + (-t26 * mrSges(4,1) + t25 * mrSges(4,3) + Ifges(4,4) * t46 + Ifges(4,6) * t48 + Ifges(4,2) * t45 / 0.2e1) * t45 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t29 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + m(3) * (t41 ^ 2 + t42 ^ 2 + t47 ^ 2) / 0.2e1 + m(4) * (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t20 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + (V_base(3) * mrSges(2,2) - t61 * mrSges(2,3) + Ifges(2,5) * t70 + Ifges(2,1) * t64 / 0.2e1) * t64 + (t47 * mrSges(3,2) - t41 * mrSges(3,3) + Ifges(3,5) * t60 + Ifges(3,1) * t54 / 0.2e1) * t54 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t33 + Ifges(6,1) * t28 / 0.2e1) * t28 + (-V_base(3) * mrSges(2,1) + t62 * mrSges(2,3) + Ifges(2,4) * t64 + Ifges(2,6) * t70 + Ifges(2,2) * t63 / 0.2e1) * t63 + (-t20 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t35 + Ifges(5,6) * t44 + Ifges(5,2) * t34 / 0.2e1) * t34 + (t26 * mrSges(4,2) - t24 * mrSges(4,3) + Ifges(4,5) * t48 + Ifges(4,1) * t46 / 0.2e1) * t46 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(2) * (t61 ^ 2 + t62 ^ 2 + t87) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t87) / 0.2e1;
T  = t8;

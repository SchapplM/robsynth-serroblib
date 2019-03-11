% Calculate Gravitation load on the joints for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:17
% EndTime: 2019-03-08 21:40:19
% DurationCPUTime: 0.87s
% Computational Cost: add. (444->146), mult. (1086->217), div. (0->0), fcn. (1280->10), ass. (0->69)
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t95 = pkin(3) * t49 + qJ(4) * t46;
t94 = rSges(6,3) + pkin(9);
t73 = rSges(7,3) + qJ(6) + pkin(9);
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t68 = cos(pkin(10));
t69 = cos(pkin(6));
t55 = t69 * t68;
t67 = sin(pkin(10));
t27 = t47 * t67 - t50 * t55;
t54 = t69 * t67;
t29 = t47 * t68 + t50 * t54;
t93 = g(1) * t29 + g(2) * t27;
t28 = t47 * t55 + t50 * t67;
t43 = sin(pkin(6));
t62 = t43 * t68;
t13 = t28 * t49 - t46 * t62;
t30 = -t47 * t54 + t50 * t68;
t61 = t43 * t67;
t15 = t30 * t49 + t46 * t61;
t79 = t43 * t47;
t32 = t46 * t69 + t49 * t79;
t92 = g(1) * t15 + g(2) * t13 + g(3) * t32;
t45 = sin(qJ(5));
t77 = t45 * t46;
t91 = pkin(5) * t77 + t49 * t73;
t90 = m(7) * t92;
t83 = g(3) * t43;
t82 = rSges(5,1) + pkin(8);
t81 = rSges(4,3) + pkin(8);
t78 = t43 * t50;
t76 = t45 * t50;
t48 = cos(qJ(5));
t75 = t46 * t48;
t74 = t48 * t50;
t72 = pkin(2) * t78 + pkin(8) * t79;
t70 = rSges(5,3) + qJ(4);
t66 = -m(5) - m(6) - m(7);
t24 = t27 * pkin(2);
t64 = -t95 * t27 - t24;
t25 = t29 * pkin(2);
t63 = -t95 * t29 - t25;
t60 = t95 * t78 + t72;
t59 = rSges(4,1) * t49 - rSges(4,2) * t46;
t58 = rSges(5,2) * t49 - rSges(5,3) * t46;
t12 = t28 * t46 + t49 * t62;
t2 = t12 * t48 - t27 * t45;
t14 = t30 * t46 - t49 * t61;
t4 = t14 * t48 - t29 * t45;
t57 = pkin(8) * t28 + t64;
t56 = pkin(8) * t30 + t63;
t31 = t46 * t79 - t49 * t69;
t16 = t31 * t48 + t43 * t76;
t42 = pkin(5) * t48 + pkin(4);
t26 = t31 * pkin(3);
t19 = (t46 * t76 + t47 * t48) * t43;
t18 = (-t45 * t47 + t46 * t74) * t43;
t17 = -t31 * t45 + t43 * t74;
t11 = t14 * pkin(3);
t10 = t12 * pkin(3);
t9 = -t29 * t77 + t30 * t48;
t8 = -t29 * t75 - t30 * t45;
t7 = -t27 * t77 + t28 * t48;
t6 = -t27 * t75 - t28 * t45;
t5 = -t14 * t45 - t29 * t48;
t3 = -t12 * t45 - t27 * t48;
t1 = [(-m(2) - m(3) - m(4) + t66) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t29 - rSges(3,2) * t30) + g(2) * (-rSges(3,1) * t27 - rSges(3,2) * t28) + (rSges(3,1) * t50 - rSges(3,2) * t47) * t83) - m(4) * (g(1) * (-t29 * t59 + t30 * t81 - t25) + g(2) * (-t27 * t59 + t28 * t81 - t24) + g(3) * t72 + (rSges(4,3) * t47 + t50 * t59) * t83) - m(5) * (g(1) * (t29 * t58 + t30 * t82 + t63) + g(2) * (t27 * t58 + t28 * t82 + t64) + g(3) * t60 + (rSges(5,1) * t47 - t50 * t58) * t83) - m(6) * (g(1) * (rSges(6,1) * t9 + rSges(6,2) * t8 + pkin(4) * t30 + t56) + g(2) * (rSges(6,1) * t7 + rSges(6,2) * t6 + pkin(4) * t28 + t57) + g(3) * (t19 * rSges(6,1) + t18 * rSges(6,2) + pkin(4) * t79 + t60) + (g(3) * t78 - t93) * t49 * t94) - m(7) * (g(1) * (rSges(7,1) * t9 + rSges(7,2) * t8 + t30 * t42 + t56) + g(2) * (rSges(7,1) * t7 + rSges(7,2) * t6 + t28 * t42 + t57) + g(3) * (t19 * rSges(7,1) + t18 * rSges(7,2) + t60) + (t42 * t47 + t91 * t50) * t83 - t93 * t91) -m(4) * (g(1) * (-rSges(4,1) * t14 - rSges(4,2) * t15) + g(2) * (-rSges(4,1) * t12 - rSges(4,2) * t13) + g(3) * (-rSges(4,1) * t31 - rSges(4,2) * t32)) - m(5) * (g(1) * (rSges(5,2) * t14 + t15 * t70 - t11) + g(2) * (rSges(5,2) * t12 + t13 * t70 - t10) + g(3) * (rSges(5,2) * t31 + t32 * t70 - t26)) - m(7) * (g(1) * (-t73 * t14 - t11) + g(2) * (-t73 * t12 - t10) + g(3) * (-t73 * t31 - t26)) - (rSges(7,2) * t48 + qJ(4) + (rSges(7,1) + pkin(5)) * t45) * t90 + (-g(1) * (-t94 * t14 - t11) - g(2) * (-t94 * t12 - t10) - g(3) * (-t94 * t31 - t26) - t92 * (rSges(6,1) * t45 + rSges(6,2) * t48 + qJ(4))) * m(6), t66 * (g(1) * t14 + g(2) * t12 + g(3) * t31) -m(6) * (g(1) * (rSges(6,1) * t4 + rSges(6,2) * t5) + g(2) * (rSges(6,1) * t2 + rSges(6,2) * t3) + g(3) * (rSges(6,1) * t16 + rSges(6,2) * t17)) + (-g(1) * (rSges(7,1) * t4 + rSges(7,2) * t5) - g(2) * (rSges(7,1) * t2 + rSges(7,2) * t3) - g(3) * (t16 * rSges(7,1) + t17 * rSges(7,2)) - (g(1) * t4 + g(2) * t2 + g(3) * t16) * pkin(5)) * m(7), -t90];
taug  = t1(:);

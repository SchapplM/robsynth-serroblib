% Calculate potential energy for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:34
% EndTime: 2019-03-09 07:27:34
% DurationCPUTime: 0.56s
% Computational Cost: add. (547->143), mult. (1191->191), div. (0->0), fcn. (1486->16), ass. (0->64)
t50 = sin(pkin(13));
t53 = cos(pkin(13));
t62 = cos(qJ(1));
t55 = cos(pkin(6));
t59 = sin(qJ(1));
t81 = t55 * t59;
t34 = -t50 * t62 - t53 * t81;
t51 = sin(pkin(7));
t54 = cos(pkin(7));
t52 = sin(pkin(6));
t84 = t52 * t59;
t24 = -t34 * t51 + t54 * t84;
t85 = t52 * t53;
t31 = -t51 * t85 + t54 * t55;
t93 = pkin(10) + rSges(5,3);
t92 = pkin(12) + rSges(7,3);
t91 = cos(qJ(3));
t80 = t55 * t62;
t32 = -t50 * t59 + t53 * t80;
t83 = t52 * t62;
t23 = -t32 * t51 - t54 * t83;
t57 = sin(qJ(4));
t90 = t23 * t57;
t89 = t24 * t57;
t88 = t31 * t57;
t86 = t50 * t52;
t79 = qJ(2) * t52;
t78 = pkin(8) + r_base(3);
t77 = t59 * pkin(1) + r_base(2);
t74 = t51 * t91;
t73 = t54 * t91;
t72 = t55 * qJ(2) + t78;
t71 = t62 * pkin(1) + t59 * t79 + r_base(1);
t70 = t52 * t74;
t35 = -t50 * t81 + t53 * t62;
t69 = t35 * pkin(2) + t24 * pkin(9) + t71;
t68 = pkin(2) * t86 + t31 * pkin(9) + t72;
t58 = sin(qJ(3));
t13 = -t34 * t73 + t35 * t58 - t59 * t70;
t14 = t35 * t91 + (t34 * t54 + t51 * t84) * t58;
t61 = cos(qJ(4));
t43 = pkin(4) * t61 + pkin(3);
t63 = -pkin(11) - pkin(10);
t67 = pkin(4) * t89 - t13 * t63 + t14 * t43 + t69;
t21 = -t55 * t74 + t58 * t86 - t73 * t85;
t22 = t55 * t51 * t58 + (t53 * t54 * t58 + t91 * t50) * t52;
t66 = pkin(4) * t88 - t21 * t63 + t22 * t43 + t68;
t33 = t50 * t80 + t53 * t59;
t65 = t33 * pkin(2) + t23 * pkin(9) - t62 * t79 + t77;
t11 = -t32 * t73 + t33 * t58 + t62 * t70;
t12 = t33 * t91 + (t32 * t54 - t51 * t83) * t58;
t64 = pkin(4) * t90 - t11 * t63 + t12 * t43 + t65;
t60 = cos(qJ(6));
t56 = sin(qJ(6));
t49 = qJ(4) + qJ(5);
t46 = cos(t49);
t45 = sin(t49);
t8 = t22 * t46 + t31 * t45;
t7 = t22 * t45 - t31 * t46;
t4 = t14 * t46 + t24 * t45;
t3 = t14 * t45 - t24 * t46;
t2 = t12 * t46 + t23 * t45;
t1 = t12 * t45 - t23 * t46;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t62 - rSges(2,2) * t59 + r_base(1)) + g(2) * (rSges(2,1) * t59 + rSges(2,2) * t62 + r_base(2)) + g(3) * (rSges(2,3) + t78)) - m(3) * (g(1) * (rSges(3,1) * t35 + rSges(3,2) * t34 + t71) + g(2) * (t33 * rSges(3,1) + t32 * rSges(3,2) + t77) + g(3) * (t55 * rSges(3,3) + t72) + (g(1) * rSges(3,3) * t59 + g(3) * (rSges(3,1) * t50 + rSges(3,2) * t53) + g(2) * (-rSges(3,3) - qJ(2)) * t62) * t52) - m(4) * (g(1) * (rSges(4,1) * t14 - rSges(4,2) * t13 + rSges(4,3) * t24 + t69) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t23 * rSges(4,3) + t65) + g(3) * (rSges(4,1) * t22 - rSges(4,2) * t21 + rSges(4,3) * t31 + t68)) - m(5) * (g(1) * (t14 * pkin(3) + (t14 * t61 + t89) * rSges(5,1) + (-t14 * t57 + t24 * t61) * rSges(5,2) + t93 * t13 + t69) + g(2) * (t12 * pkin(3) + (t12 * t61 + t90) * rSges(5,1) + (-t12 * t57 + t23 * t61) * rSges(5,2) + t93 * t11 + t65) + g(3) * (t22 * pkin(3) + (t22 * t61 + t88) * rSges(5,1) + (-t22 * t57 + t31 * t61) * rSges(5,2) + t93 * t21 + t68)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t13 + t67) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t11 + t64) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + rSges(6,3) * t21 + t66)) - m(7) * (g(1) * (t4 * pkin(5) + (t13 * t56 + t4 * t60) * rSges(7,1) + (t13 * t60 - t4 * t56) * rSges(7,2) + t92 * t3 + t67) + g(2) * (t2 * pkin(5) + (t11 * t56 + t2 * t60) * rSges(7,1) + (t11 * t60 - t2 * t56) * rSges(7,2) + t92 * t1 + t64) + g(3) * (t8 * pkin(5) + (t21 * t56 + t60 * t8) * rSges(7,1) + (t21 * t60 - t56 * t8) * rSges(7,2) + t92 * t7 + t66));
U  = t5;

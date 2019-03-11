% Calculate potential energy for
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:37
% EndTime: 2019-03-08 18:59:37
% DurationCPUTime: 0.58s
% Computational Cost: add. (547->143), mult. (1191->191), div. (0->0), fcn. (1486->16), ass. (0->64)
t52 = sin(pkin(7));
t56 = cos(pkin(7));
t57 = cos(pkin(6));
t53 = sin(pkin(6));
t54 = cos(pkin(13));
t83 = t53 * t54;
t31 = -t52 * t83 + t56 * t57;
t50 = sin(pkin(13));
t55 = cos(pkin(12));
t51 = sin(pkin(12));
t85 = t51 * t57;
t34 = -t50 * t55 - t54 * t85;
t82 = t53 * t56;
t24 = -t34 * t52 + t51 * t82;
t93 = pkin(9) + rSges(5,3);
t92 = pkin(11) + rSges(7,3);
t91 = cos(qJ(3));
t81 = t55 * t57;
t32 = -t50 * t51 + t54 * t81;
t23 = -t32 * t52 - t55 * t82;
t59 = sin(qJ(4));
t90 = t23 * t59;
t89 = t24 * t59;
t88 = t31 * t59;
t86 = t50 * t53;
t84 = t52 * t53;
t79 = qJ(2) * t53;
t78 = t51 * pkin(1) + r_base(2);
t75 = qJ(1) + r_base(3);
t74 = t52 * t91;
t73 = t56 * t91;
t72 = t55 * pkin(1) + t51 * t79 + r_base(1);
t71 = t57 * qJ(2) + t75;
t70 = t53 * t74;
t35 = -t50 * t85 + t54 * t55;
t69 = t35 * pkin(2) + t24 * pkin(8) + t72;
t68 = pkin(2) * t86 + t31 * pkin(8) + t71;
t60 = sin(qJ(3));
t13 = -t34 * t73 + t35 * t60 - t51 * t70;
t14 = t35 * t91 + (t34 * t56 + t51 * t84) * t60;
t62 = cos(qJ(4));
t43 = pkin(4) * t62 + pkin(3);
t63 = -pkin(10) - pkin(9);
t67 = pkin(4) * t89 - t13 * t63 + t14 * t43 + t69;
t21 = -t57 * t74 + t60 * t86 - t73 * t83;
t22 = t57 * t52 * t60 + (t54 * t56 * t60 + t91 * t50) * t53;
t66 = pkin(4) * t88 - t21 * t63 + t22 * t43 + t68;
t33 = t50 * t81 + t51 * t54;
t65 = t33 * pkin(2) + t23 * pkin(8) - t55 * t79 + t78;
t11 = -t32 * t73 + t33 * t60 + t55 * t70;
t12 = t33 * t91 + (t32 * t56 - t55 * t84) * t60;
t64 = pkin(4) * t90 - t11 * t63 + t12 * t43 + t65;
t61 = cos(qJ(6));
t58 = sin(qJ(6));
t49 = qJ(4) + qJ(5);
t46 = cos(t49);
t45 = sin(t49);
t8 = t22 * t46 + t31 * t45;
t7 = t22 * t45 - t31 * t46;
t4 = t14 * t46 + t24 * t45;
t3 = t14 * t45 - t24 * t46;
t2 = t12 * t46 + t23 * t45;
t1 = t12 * t45 - t23 * t46;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t55 - rSges(2,2) * t51 + r_base(1)) + g(2) * (rSges(2,1) * t51 + rSges(2,2) * t55 + r_base(2)) + g(3) * (rSges(2,3) + t75)) - m(3) * (g(1) * (rSges(3,1) * t35 + rSges(3,2) * t34 + t72) + g(2) * (t33 * rSges(3,1) + t32 * rSges(3,2) + t78) + g(3) * (t57 * rSges(3,3) + t71) + (g(1) * rSges(3,3) * t51 + g(3) * (rSges(3,1) * t50 + rSges(3,2) * t54) + g(2) * (-rSges(3,3) - qJ(2)) * t55) * t53) - m(4) * (g(1) * (rSges(4,1) * t14 - rSges(4,2) * t13 + rSges(4,3) * t24 + t69) + g(2) * (t12 * rSges(4,1) - t11 * rSges(4,2) + t23 * rSges(4,3) + t65) + g(3) * (rSges(4,1) * t22 - rSges(4,2) * t21 + rSges(4,3) * t31 + t68)) - m(5) * (g(1) * (t14 * pkin(3) + (t14 * t62 + t89) * rSges(5,1) + (-t14 * t59 + t24 * t62) * rSges(5,2) + t93 * t13 + t69) + g(2) * (t12 * pkin(3) + (t12 * t62 + t90) * rSges(5,1) + (-t12 * t59 + t23 * t62) * rSges(5,2) + t93 * t11 + t65) + g(3) * (t22 * pkin(3) + (t22 * t62 + t88) * rSges(5,1) + (-t22 * t59 + t31 * t62) * rSges(5,2) + t93 * t21 + t68)) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t13 + t67) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + rSges(6,3) * t11 + t64) + g(3) * (rSges(6,1) * t8 - rSges(6,2) * t7 + rSges(6,3) * t21 + t66)) - m(7) * (g(1) * (t4 * pkin(5) + (t13 * t58 + t4 * t61) * rSges(7,1) + (t13 * t61 - t4 * t58) * rSges(7,2) + t92 * t3 + t67) + g(2) * (t2 * pkin(5) + (t11 * t58 + t2 * t61) * rSges(7,1) + (t11 * t61 - t2 * t58) * rSges(7,2) + t92 * t1 + t64) + g(3) * (t8 * pkin(5) + (t21 * t58 + t61 * t8) * rSges(7,1) + (t21 * t61 - t58 * t8) * rSges(7,2) + t92 * t7 + t66));
U  = t5;

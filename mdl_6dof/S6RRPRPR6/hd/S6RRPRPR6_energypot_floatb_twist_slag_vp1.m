% Calculate potential energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:01
% EndTime: 2019-03-09 10:38:01
% DurationCPUTime: 0.50s
% Computational Cost: add. (393->119), mult. (818->142), div. (0->0), fcn. (1011->12), ass. (0->55)
t35 = sin(pkin(11));
t40 = sin(qJ(2));
t44 = cos(qJ(2));
t61 = cos(pkin(11));
t24 = -t40 * t35 + t44 * t61;
t74 = pkin(2) * t40;
t73 = rSges(6,1) + pkin(9);
t71 = rSges(5,3) + pkin(9);
t70 = pkin(10) + rSges(7,3);
t36 = sin(pkin(6));
t41 = sin(qJ(1));
t68 = t41 * t36;
t67 = t41 * t40;
t66 = t41 * t44;
t45 = cos(qJ(1));
t65 = t45 * t36;
t64 = t45 * t40;
t63 = t45 * t44;
t62 = rSges(6,3) + qJ(5);
t60 = pkin(7) + r_base(3);
t32 = t44 * pkin(2) + pkin(1);
t59 = t45 * t32 + r_base(1);
t37 = cos(pkin(6));
t58 = t37 * pkin(8) + t60;
t22 = t37 * t74 + (-pkin(8) - qJ(3)) * t36;
t56 = t45 * t22 + t41 * t32 + r_base(2);
t23 = -t44 * t35 - t40 * t61;
t21 = t23 * t37;
t10 = -t45 * t21 + t41 * t24;
t55 = t10 * pkin(3) + t56;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t4 = t10 * t43 - t39 * t65;
t54 = t4 * pkin(4) + t55;
t53 = t37 * qJ(3) + t36 * t74 + t58;
t12 = t41 * t21 + t45 * t24;
t52 = t12 * pkin(3) - t41 * t22 + t59;
t38 = sin(qJ(6));
t42 = cos(qJ(6));
t51 = t38 * rSges(7,1) + t42 * rSges(7,2) + qJ(5);
t20 = t23 * t36;
t50 = -t20 * pkin(3) + t53;
t49 = t42 * rSges(7,1) - t38 * rSges(7,2) + pkin(5) + pkin(9);
t6 = t12 * t43 + t39 * t68;
t48 = t6 * pkin(4) + t52;
t15 = -t20 * t43 + t37 * t39;
t47 = t15 * pkin(4) + t50;
t46 = t24 * t37;
t19 = t24 * t36;
t14 = -t20 * t39 - t37 * t43;
t11 = t45 * t23 - t41 * t46;
t9 = t41 * t23 + t45 * t46;
t5 = t12 * t39 - t43 * t68;
t3 = t10 * t39 + t43 * t65;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t45 * rSges(2,1) - t41 * rSges(2,2) + r_base(1)) + g(2) * (t41 * rSges(2,1) + t45 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t60)) - m(3) * (g(1) * (t45 * pkin(1) + r_base(1) + (-t37 * t67 + t63) * rSges(3,1) + (-t37 * t66 - t64) * rSges(3,2)) + g(2) * (t41 * pkin(1) + r_base(2) + (t37 * t64 + t66) * rSges(3,1) + (t37 * t63 - t67) * rSges(3,2)) + g(3) * (t37 * rSges(3,3) + t58) + (g(3) * (rSges(3,1) * t40 + rSges(3,2) * t44) + (g(1) * t41 - g(2) * t45) * (rSges(3,3) + pkin(8))) * t36) - m(4) * (g(1) * (t12 * rSges(4,1) + t11 * rSges(4,2) + (rSges(4,3) * t36 - t22) * t41 + t59) + g(2) * (t10 * rSges(4,1) + t9 * rSges(4,2) - rSges(4,3) * t65 + t56) + g(3) * (-t20 * rSges(4,1) + t19 * rSges(4,2) + t37 * rSges(4,3) + t53)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) - t11 * t71 + t52) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) - t71 * t9 + t55) + g(3) * (t15 * rSges(5,1) - t14 * rSges(5,2) - t19 * t71 + t50)) - m(6) * (g(1) * (-t6 * rSges(6,2) - t11 * t73 + t5 * t62 + t48) + g(2) * (-t4 * rSges(6,2) + t3 * t62 - t73 * t9 + t54) + g(3) * (-t15 * rSges(6,2) + t14 * t62 - t19 * t73 + t47)) - m(7) * (g(1) * (-t11 * t49 + t5 * t51 + t6 * t70 + t48) + g(2) * (t3 * t51 + t4 * t70 - t49 * t9 + t54) + g(3) * (t14 * t51 + t15 * t70 - t19 * t49 + t47));
U  = t1;

% Calculate Gravitation load on the joints for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:47
% EndTime: 2019-03-09 07:04:49
% DurationCPUTime: 0.55s
% Computational Cost: add. (553->111), mult. (387->139), div. (0->0), fcn. (332->12), ass. (0->65)
t93 = rSges(7,3) + pkin(10);
t38 = pkin(11) + qJ(3);
t34 = qJ(4) + t38;
t31 = qJ(5) + t34;
t25 = sin(t31);
t26 = cos(t31);
t42 = sin(qJ(6));
t82 = rSges(7,2) * t42;
t92 = t25 * t82 + t26 * t93;
t90 = t26 * rSges(6,1) - t25 * rSges(6,2);
t28 = sin(t34);
t29 = cos(t34);
t65 = t29 * rSges(5,1) - t28 * rSges(5,2);
t43 = sin(qJ(1));
t45 = cos(qJ(1));
t89 = g(1) * t45 + g(2) * t43;
t52 = t26 * pkin(5) + t93 * t25;
t44 = cos(qJ(6));
t83 = rSges(7,1) * t44;
t88 = (-pkin(5) - t83) * t25;
t32 = sin(t38);
t87 = pkin(3) * t32;
t86 = pkin(4) * t28;
t40 = cos(pkin(11));
t30 = t40 * pkin(2) + pkin(1);
t77 = t43 * t42;
t76 = t43 * t44;
t75 = t45 * t42;
t74 = t45 * t44;
t41 = -pkin(7) - qJ(2);
t73 = rSges(4,3) - t41;
t37 = -pkin(8) + t41;
t72 = rSges(5,3) - t37;
t35 = -pkin(9) + t37;
t71 = rSges(6,3) - t35;
t70 = rSges(3,3) + qJ(2);
t33 = cos(t38);
t27 = pkin(3) * t33;
t12 = t27 + t30;
t68 = t92 * t43;
t67 = t92 * t45;
t24 = pkin(4) * t29;
t63 = t24 + t90;
t62 = t33 * rSges(4,1) - t32 * rSges(4,2);
t60 = -rSges(5,1) * t28 - rSges(5,2) * t29;
t58 = -rSges(6,1) * t25 - rSges(6,2) * t26;
t57 = rSges(3,1) * t40 - rSges(3,2) * sin(pkin(11)) + pkin(1);
t56 = t30 + t62;
t55 = t12 + t65;
t54 = t58 * t43;
t53 = t58 * t45;
t51 = t52 + (-t82 + t83) * t26;
t49 = t24 + t51;
t48 = g(1) * t67 + g(2) * t68;
t47 = t89 * t88;
t9 = -t86 - t87;
t8 = t45 * t9;
t7 = t43 * t9;
t6 = t26 * t74 + t77;
t5 = -t26 * t75 + t76;
t4 = -t26 * t76 + t75;
t3 = t26 * t77 + t74;
t2 = t24 + t12;
t1 = t45 * t2;
t10 = [-m(2) * (g(1) * (-t43 * rSges(2,1) - t45 * rSges(2,2)) + g(2) * (t45 * rSges(2,1) - t43 * rSges(2,2))) - m(3) * ((g(1) * t70 + g(2) * t57) * t45 + (-g(1) * t57 + g(2) * t70) * t43) - m(4) * ((g(1) * t73 + g(2) * t56) * t45 + (-g(1) * t56 + g(2) * t73) * t43) - m(5) * ((g(1) * t72 + g(2) * t55) * t45 + (-g(1) * t55 + g(2) * t72) * t43) - m(6) * (g(2) * t1 + (g(1) * t71 + g(2) * t90) * t45 + (g(1) * (-t2 - t90) + g(2) * t71) * t43) - m(7) * (g(1) * (t4 * rSges(7,1) + t3 * rSges(7,2)) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t1) + (-g(1) * t35 + g(2) * t52) * t45 + (g(1) * (-t2 - t52) - g(2) * t35) * t43) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t43 - g(2) * t45) -m(4) * (g(3) * t62 + t89 * (-rSges(4,1) * t32 - rSges(4,2) * t33)) - m(5) * (g(3) * (t27 + t65) + t89 * (t60 - t87)) - m(6) * (g(1) * (t8 + t53) + g(2) * (t7 + t54) + g(3) * (t27 + t63)) - m(7) * (g(1) * (t8 + t67) + g(2) * (t7 + t68) + g(3) * (t27 + t49) + t47) -m(7) * t48 + (-m(5) * t65 - m(6) * t63 - m(7) * t49) * g(3) + t89 * (-m(5) * t60 - m(6) * (t58 - t86) - m(7) * (-t86 + t88)) -m(6) * (g(1) * t53 + g(2) * t54 + g(3) * t90) - m(7) * (g(3) * t51 + t47 + t48) -m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * (-t3 * rSges(7,1) + t4 * rSges(7,2)) + g(3) * (-rSges(7,1) * t42 - rSges(7,2) * t44) * t25)];
taug  = t10(:);

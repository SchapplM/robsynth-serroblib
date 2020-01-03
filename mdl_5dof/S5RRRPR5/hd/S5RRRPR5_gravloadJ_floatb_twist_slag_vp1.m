% Calculate Gravitation load on the joints for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:01
% EndTime: 2019-12-31 21:13:03
% DurationCPUTime: 0.47s
% Computational Cost: add. (319->93), mult. (304->123), div. (0->0), fcn. (263->10), ass. (0->54)
t78 = rSges(6,3) + pkin(8);
t32 = qJ(2) + qJ(3);
t27 = pkin(9) + t32;
t23 = sin(t27);
t24 = cos(t27);
t33 = sin(qJ(5));
t67 = rSges(6,2) * t33;
t77 = t23 * t67 + t24 * t78;
t47 = t24 * rSges(5,1) - t23 * rSges(5,2);
t28 = sin(t32);
t29 = cos(t32);
t52 = t29 * rSges(4,1) - t28 * rSges(4,2);
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t75 = g(1) * t38 + g(2) * t35;
t43 = t24 * pkin(4) + t78 * t23;
t36 = cos(qJ(5));
t68 = rSges(6,1) * t36;
t74 = (-pkin(4) - t68) * t23;
t39 = -pkin(7) - pkin(6);
t73 = pkin(3) * t28;
t34 = sin(qJ(2));
t70 = t34 * pkin(2);
t69 = rSges(3,3) + pkin(6);
t37 = cos(qJ(2));
t30 = t37 * pkin(2);
t26 = t30 + pkin(1);
t62 = t35 * t33;
t61 = t35 * t36;
t60 = t38 * t33;
t59 = t38 * t36;
t58 = rSges(4,3) - t39;
t31 = -qJ(4) + t39;
t57 = rSges(5,3) - t31;
t56 = t77 * t35;
t55 = t77 * t38;
t25 = pkin(3) * t29;
t51 = t25 + t47;
t50 = t37 * rSges(3,1) - t34 * rSges(3,2);
t48 = -rSges(4,1) * t28 - rSges(4,2) * t29;
t46 = -rSges(5,1) * t23 - rSges(5,2) * t24;
t45 = pkin(1) + t50;
t44 = t26 + t52;
t41 = t25 + t43 + (-t67 + t68) * t24;
t11 = -t70 - t73;
t10 = t25 + t26;
t7 = t38 * t11;
t6 = t35 * t11;
t5 = t38 * t10;
t4 = t24 * t59 + t62;
t3 = -t24 * t60 + t61;
t2 = -t24 * t61 + t60;
t1 = t24 * t62 + t59;
t8 = [-m(2) * (g(1) * (-t35 * rSges(2,1) - t38 * rSges(2,2)) + g(2) * (t38 * rSges(2,1) - t35 * rSges(2,2))) - m(3) * ((g(1) * t69 + g(2) * t45) * t38 + (-g(1) * t45 + g(2) * t69) * t35) - m(4) * ((g(1) * t58 + g(2) * t44) * t38 + (-g(1) * t44 + g(2) * t58) * t35) - m(5) * (g(2) * t5 + (g(1) * t57 + g(2) * t47) * t38 + (g(1) * (-t10 - t47) + g(2) * t57) * t35) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (-g(1) * t31 + g(2) * t43) * t38 + (g(1) * (-t10 - t43) - g(2) * t31) * t35), -m(3) * (g(3) * t50 + t75 * (-rSges(3,1) * t34 - rSges(3,2) * t37)) - m(4) * (g(3) * (t30 + t52) + t75 * (t48 - t70)) - m(5) * (g(1) * (t46 * t38 + t7) + g(2) * (t46 * t35 + t6) + g(3) * (t30 + t51)) - m(6) * (g(1) * (t7 + t55) + g(2) * (t6 + t56) + g(3) * (t30 + t41) + t75 * t74), -m(6) * (g(1) * t55 + g(2) * t56) + (-m(4) * t52 - m(5) * t51 - m(6) * t41) * g(3) + t75 * (-m(4) * t48 - m(5) * (t46 - t73) - m(6) * (-t73 + t74)), (-m(5) - m(6)) * (g(1) * t35 - g(2) * t38), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t33 - rSges(6,2) * t36) * t23)];
taug = t8(:);

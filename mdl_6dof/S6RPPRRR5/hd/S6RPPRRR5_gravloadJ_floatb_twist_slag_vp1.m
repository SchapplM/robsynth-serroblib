% Calculate Gravitation load on the joints for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:25
% EndTime: 2019-03-09 02:28:26
% DurationCPUTime: 0.45s
% Computational Cost: add. (219->97), mult. (293->130), div. (0->0), fcn. (250->8), ass. (0->51)
t32 = sin(qJ(6));
t35 = cos(qJ(6));
t80 = rSges(7,1) * t35 - rSges(7,2) * t32;
t79 = rSges(7,3) + pkin(9);
t31 = qJ(4) + qJ(5);
t26 = cos(t31);
t78 = t80 * t26;
t77 = t79 * t26;
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t76 = g(1) * t37 + g(2) * t34;
t25 = sin(t31);
t75 = t25 * pkin(5) - t77;
t36 = cos(qJ(4));
t74 = pkin(4) * t36;
t73 = g(1) * t34;
t70 = g(2) * t37;
t33 = sin(qJ(4));
t68 = t33 * pkin(4);
t67 = -rSges(5,3) - pkin(7);
t63 = t25 * t34;
t62 = t25 * t37;
t61 = t26 * t34;
t60 = t26 * t37;
t59 = t34 * t32;
t58 = t34 * t35;
t57 = t37 * t32;
t56 = t37 * t35;
t55 = -pkin(1) - qJ(3);
t29 = t37 * qJ(2);
t38 = -pkin(8) - pkin(7);
t54 = t37 * t38 + t29;
t53 = t37 * pkin(1) + t34 * qJ(2);
t50 = t37 * qJ(3) + t53;
t49 = -m(4) - m(5) - m(6) - m(7);
t48 = rSges(6,1) * t61 - rSges(6,2) * t63;
t47 = rSges(6,1) * t60 - rSges(6,2) * t62;
t46 = t34 * t38 + t37 * t68 + t50;
t44 = t33 * rSges(5,1) + t36 * rSges(5,2);
t43 = -t25 * rSges(6,1) - t26 * rSges(6,2);
t42 = t43 - t68;
t41 = pkin(5) * t61 + t78 * t34 + t79 * t63;
t40 = pkin(5) * t60 + t78 * t37 + t79 * t62;
t39 = t77 + (-pkin(5) - t80) * t25;
t20 = t37 * t74;
t18 = t34 * t74;
t4 = t25 * t56 - t59;
t3 = -t25 * t57 - t58;
t2 = -t25 * t58 - t57;
t1 = t25 * t59 - t56;
t5 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - t37 * rSges(2,2)) + g(2) * (t37 * rSges(2,1) - t34 * rSges(2,2))) - m(3) * (g(1) * (t37 * rSges(3,3) + t29 + (rSges(3,2) - pkin(1)) * t34) + g(2) * (-t37 * rSges(3,2) + t34 * rSges(3,3) + t53)) - m(4) * (g(1) * (t37 * rSges(4,2) + t29) + g(2) * (t37 * rSges(4,3) + t50) + (g(1) * (-rSges(4,3) + t55) + g(2) * rSges(4,2)) * t34) - m(5) * (g(1) * t29 + g(2) * t50 + (g(1) * t67 + g(2) * t44) * t37 + (g(1) * (-t44 + t55) + g(2) * t67) * t34) - m(6) * (g(1) * (-t37 * rSges(6,3) + t54) + g(2) * (rSges(6,1) * t62 + rSges(6,2) * t60 + t46) + (g(1) * (t42 + t55) - g(2) * rSges(6,3)) * t34) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t54) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t46) + t75 * t70 + (t55 - t68 - t75) * t73) (-m(3) + t49) * (-t70 + t73) t49 * t76, -m(5) * (-g(3) * t44 + t76 * (rSges(5,1) * t36 - rSges(5,2) * t33)) - m(6) * (g(1) * (t20 + t47) + g(2) * (t18 + t48) + g(3) * t42) - m(7) * (g(1) * (t20 + t40) + g(2) * (t18 + t41) + g(3) * (t39 - t68)) -m(6) * (g(1) * t47 + g(2) * t48 + g(3) * t43) - m(7) * (g(1) * t40 + g(2) * t41 + g(3) * t39) -m(7) * (g(1) * (t3 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (-t1 * rSges(7,1) + t2 * rSges(7,2)) + g(3) * (-rSges(7,1) * t32 - rSges(7,2) * t35) * t26)];
taug  = t5(:);

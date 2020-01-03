% Calculate Gravitation load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:05
% EndTime: 2020-01-03 12:11:06
% DurationCPUTime: 0.30s
% Computational Cost: add. (290->84), mult. (234->108), div. (0->0), fcn. (186->8), ass. (0->47)
t74 = rSges(4,3) + pkin(7);
t46 = -pkin(8) - pkin(7);
t73 = rSges(5,3) - t46;
t72 = rSges(6,3) + qJ(5) - t46;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t71 = t44 * rSges(4,1) - t42 * rSges(4,2);
t70 = pkin(2) + t71;
t40 = qJ(3) + qJ(4);
t32 = sin(t40);
t69 = pkin(4) * t32;
t41 = qJ(1) + qJ(2);
t33 = sin(t41);
t68 = g(2) * t33;
t67 = t42 * pkin(3);
t34 = cos(t40);
t35 = cos(t41);
t63 = t34 * t35;
t64 = t32 * t35;
t66 = rSges(6,1) * t64 + rSges(6,2) * t63;
t37 = t44 * pkin(3);
t31 = t37 + pkin(2);
t65 = t32 * t33;
t24 = t34 * rSges(5,1);
t23 = t34 * rSges(6,1);
t60 = rSges(5,1) * t64 + rSges(5,2) * t63;
t59 = t33 * rSges(3,1) + t35 * rSges(3,2);
t58 = t35 * rSges(3,1) - t33 * rSges(3,2);
t57 = -t32 * rSges(5,2) + t24;
t30 = pkin(4) * t34;
t56 = -t32 * rSges(6,2) + t23 + t30;
t55 = rSges(4,1) * t42 + rSges(4,2) * t44;
t54 = -rSges(5,1) * t32 - rSges(5,2) * t34;
t53 = -rSges(6,1) * t32 - rSges(6,2) * t34;
t52 = t74 * t33 + t70 * t35;
t51 = -rSges(5,2) * t65 - t73 * t35 + (t24 + t31) * t33;
t3 = t30 + t31;
t50 = -rSges(6,2) * t65 - t72 * t35 + (t23 + t3) * t33;
t49 = rSges(5,1) * t63 - rSges(5,2) * t64 + t35 * t31 + t73 * t33;
t48 = rSges(6,1) * t63 - rSges(6,2) * t64 + t35 * t3 + t72 * t33;
t47 = t70 * t33 - t74 * t35;
t45 = cos(qJ(1));
t43 = sin(qJ(1));
t38 = t45 * pkin(1);
t36 = t43 * pkin(1);
t4 = -t67 - t69;
t1 = [-m(2) * (g(2) * (t45 * rSges(2,1) - t43 * rSges(2,2)) + g(3) * (t43 * rSges(2,1) + t45 * rSges(2,2))) - m(3) * (g(2) * (t38 + t58) + g(3) * (t36 + t59)) - m(4) * (g(2) * (t38 + t52) + g(3) * (t36 + t47)) - m(5) * (g(2) * (t38 + t49) + g(3) * (t36 + t51)) - m(6) * (g(2) * (t38 + t48) + g(3) * (t36 + t50)), -m(3) * (g(2) * t58 + g(3) * t59) - m(4) * (g(2) * t52 + g(3) * t47) - m(5) * (g(2) * t49 + g(3) * t51) - m(6) * (g(2) * t48 + g(3) * t50), -m(4) * (g(3) * t55 * t35 + g(1) * t71) - m(5) * (g(1) * (t37 + t57) + g(3) * (t35 * t67 + t60)) - m(6) * (g(1) * (t37 + t56) + g(3) * (-t35 * t4 + t66)) + (m(4) * t55 - m(5) * (t54 - t67) - m(6) * (t4 + t53)) * t68, -m(5) * (g(1) * t57 + g(3) * t60) - m(6) * (g(1) * t56 + g(3) * (pkin(4) * t64 + t66)) + (-m(5) * t54 - m(6) * (t53 - t69)) * t68, -m(6) * (-g(2) * t35 - g(3) * t33)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:40
% EndTime: 2019-12-31 21:50:41
% DurationCPUTime: 0.37s
% Computational Cost: add. (324->79), mult. (261->98), div. (0->0), fcn. (213->8), ass. (0->43)
t75 = rSges(6,1) + pkin(4);
t37 = qJ(1) + qJ(2);
t31 = sin(t37);
t76 = g(2) * t31;
t74 = rSges(4,3) + pkin(7);
t36 = qJ(3) + qJ(4);
t30 = sin(t36);
t32 = cos(t36);
t55 = t32 * rSges(5,1) - t30 * rSges(5,2);
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t73 = t40 * rSges(4,1) - t38 * rSges(4,2);
t72 = -pkin(2) - t73;
t33 = cos(t37);
t70 = g(1) * t33 + t76;
t69 = t70 * t30;
t14 = t30 * qJ(5);
t53 = t30 * rSges(6,3) + t75 * t32 + t14;
t68 = pkin(3) * t38;
t39 = sin(qJ(1));
t65 = t39 * pkin(1);
t62 = t30 * t33;
t61 = t32 * t33;
t42 = -pkin(8) - pkin(7);
t60 = t33 * t42;
t57 = qJ(5) * t32;
t56 = t33 * rSges(3,1) - t31 * rSges(3,2);
t54 = t33 * rSges(6,2) - t60;
t34 = t40 * pkin(3);
t29 = t34 + pkin(2);
t4 = t33 * t29;
t52 = t31 * rSges(6,2) + rSges(6,3) * t62 + t33 * t14 + t75 * t61 + t4;
t51 = -t31 * rSges(3,1) - t33 * rSges(3,2);
t50 = -rSges(5,1) * t30 - rSges(5,2) * t32;
t49 = g(1) * (rSges(6,3) * t61 + t33 * t57) + (rSges(6,3) * t32 + t57) * t76;
t48 = t74 * t31 - t72 * t33;
t47 = t72 * t31 + t74 * t33;
t46 = rSges(5,1) * t61 - rSges(5,2) * t62 + t4 + (rSges(5,3) - t42) * t31;
t45 = t33 * rSges(5,3) - t60 + (-t29 - t55) * t31;
t44 = (g(1) * (-t29 - t53) - g(2) * t42) * t31;
t41 = cos(qJ(1));
t35 = t41 * pkin(1);
t1 = [-m(2) * (g(1) * (-t39 * rSges(2,1) - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(1) * (t51 - t65) + g(2) * (t35 + t56)) - m(4) * (g(1) * (t47 - t65) + g(2) * (t35 + t48)) - m(5) * (g(1) * (t45 - t65) + g(2) * (t35 + t46)) - m(6) * (g(1) * (t54 - t65) + g(2) * (t35 + t52) + t44), -m(3) * (g(1) * t51 + g(2) * t56) - m(4) * (g(1) * t47 + g(2) * t48) - m(5) * (g(1) * t45 + g(2) * t46) - m(6) * (g(1) * t54 + g(2) * t52 + t44), -m(6) * t49 + (-m(4) * t73 - m(5) * (t34 + t55) - m(6) * (t34 + t53)) * g(3) + t70 * (-m(4) * (-rSges(4,1) * t38 - rSges(4,2) * t40) - m(5) * (t50 - t68) - m(6) * (-t30 * t75 - t68)), -m(5) * (g(3) * t55 + t70 * t50) - m(6) * (g(3) * t53 - t69 * t75 + t49), -m(6) * (-g(3) * t32 + t69)];
taug = t1(:);

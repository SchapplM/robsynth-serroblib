% Calculate Gravitation load on the joints for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:01
% DurationCPUTime: 0.35s
% Computational Cost: add. (303->88), mult. (324->114), div. (0->0), fcn. (299->8), ass. (0->46)
t40 = cos(qJ(3));
t37 = sin(qJ(3));
t32 = t37 * qJ(4);
t59 = t40 * pkin(3) + t32;
t78 = -t40 * pkin(4) - t59;
t77 = t40 * rSges(5,1) + t37 * rSges(5,3);
t76 = t40 * rSges(4,1) - t37 * rSges(4,2);
t75 = -rSges(6,3) - pkin(8);
t36 = sin(qJ(5));
t39 = cos(qJ(5));
t51 = t40 * t36 - t37 * t39;
t35 = qJ(1) + qJ(2);
t30 = sin(t35);
t31 = cos(t35);
t74 = g(1) * t31 + g(2) * t30;
t38 = sin(qJ(1));
t71 = t38 * pkin(1);
t69 = -rSges(5,1) - pkin(3);
t68 = t31 * t40;
t28 = t31 * pkin(7);
t61 = t31 * rSges(5,2) + t28;
t60 = t31 * pkin(2) + t30 * pkin(7);
t58 = qJ(4) * t40;
t57 = t31 * rSges(3,1) - t30 * rSges(3,2);
t56 = pkin(3) * t68 + t31 * t32 + t60;
t6 = t51 * t30;
t50 = t37 * t36 + t40 * t39;
t7 = t50 * t30;
t55 = -t6 * rSges(6,1) - t7 * rSges(6,2);
t8 = t51 * t31;
t9 = t50 * t31;
t54 = -t8 * rSges(6,1) - t9 * rSges(6,2);
t53 = -t30 * rSges(3,1) - t31 * rSges(3,2);
t52 = -rSges(6,1) * t50 + rSges(6,2) * t51;
t49 = t9 * rSges(6,1) - t8 * rSges(6,2) + pkin(4) * t68 + t56;
t48 = t30 * rSges(5,2) + t77 * t31 + t56;
t47 = t30 * rSges(4,3) + t76 * t31 + t60;
t46 = -t7 * rSges(6,1) + t6 * rSges(6,2) + t75 * t31 + t28;
t45 = t31 * rSges(4,3) + t28 + (-pkin(2) - t76) * t30;
t44 = g(1) * (-pkin(2) + t69 * t40 + (-rSges(5,3) - qJ(4)) * t37) * t30;
t43 = (g(1) * (-pkin(2) + t78) + g(2) * t75) * t30;
t41 = cos(qJ(1));
t34 = t41 * pkin(1);
t15 = t31 * t58;
t13 = t30 * t58;
t1 = [-m(2) * (g(1) * (-t38 * rSges(2,1) - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) - t38 * rSges(2,2))) - m(3) * (g(1) * (t53 - t71) + g(2) * (t34 + t57)) - m(4) * (g(1) * (t45 - t71) + g(2) * (t34 + t47)) - m(5) * (g(1) * (t61 - t71) + g(2) * (t34 + t48) + t44) - m(6) * (g(1) * (t46 - t71) + g(2) * (t34 + t49) + t43), -m(3) * (g(1) * t53 + g(2) * t57) - m(4) * (g(1) * t45 + g(2) * t47) - m(5) * (g(1) * t61 + g(2) * t48 + t44) - m(6) * (g(1) * t46 + g(2) * t49 + t43), -m(4) * g(3) * t76 - m(5) * (g(1) * t15 + g(2) * t13 + g(3) * (t59 + t77)) - m(6) * (g(1) * (t15 - t54) + g(2) * (t13 - t55) + g(3) * (-t52 - t78)) + t74 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3)) * t40 + (m(4) * rSges(4,1) - m(5) * t69 - m(6) * (-pkin(3) - pkin(4))) * t37), (-m(5) - m(6)) * (-g(3) * t40 + t74 * t37), -m(6) * (g(1) * t54 + g(2) * t55 + g(3) * t52)];
taug = t1(:);

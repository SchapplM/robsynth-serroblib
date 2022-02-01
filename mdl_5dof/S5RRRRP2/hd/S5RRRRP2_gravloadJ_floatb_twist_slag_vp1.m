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
% m [6x1]
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:48:45
% EndTime: 2022-01-20 11:48:46
% DurationCPUTime: 0.32s
% Computational Cost: add. (290->74), mult. (234->92), div. (0->0), fcn. (186->8), ass. (0->42)
t70 = rSges(4,3) + pkin(7);
t39 = -pkin(8) - pkin(7);
t69 = qJ(5) - t39 + rSges(6,3);
t68 = -t39 + rSges(5,3);
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t67 = t37 * rSges(4,1) - t35 * rSges(4,2);
t66 = -pkin(2) - t67;
t33 = qJ(3) + qJ(4);
t26 = sin(t33);
t28 = cos(t33);
t65 = -rSges(6,2) * t28 + (-rSges(6,1) - pkin(4)) * t26;
t34 = qJ(1) + qJ(2);
t27 = sin(t34);
t29 = cos(t34);
t64 = g(1) * t29 + g(2) * t27;
t60 = t35 * pkin(3);
t36 = sin(qJ(1));
t59 = t36 * pkin(1);
t30 = t37 * pkin(3);
t25 = t30 + pkin(2);
t58 = t26 * t27;
t57 = t26 * t29;
t16 = t28 * rSges(5,1);
t15 = t28 * rSges(6,1);
t56 = t28 * t29;
t53 = t29 * rSges(3,1) - t27 * rSges(3,2);
t52 = -t26 * rSges(5,2) + t16;
t24 = pkin(4) * t28;
t51 = -t26 * rSges(6,2) + t15 + t24;
t50 = -t27 * rSges(3,1) - t29 * rSges(3,2);
t49 = -rSges(5,1) * t26 - rSges(5,2) * t28;
t47 = t70 * t27 - t66 * t29;
t2 = t24 + t25;
t46 = rSges(6,2) * t58 + (-t2 - t15) * t27 + t69 * t29;
t45 = t66 * t27 + t70 * t29;
t44 = rSges(5,1) * t56 - rSges(5,2) * t57 + t29 * t25 + t68 * t27;
t43 = rSges(6,1) * t56 - rSges(6,2) * t57 + t29 * t2 + t69 * t27;
t42 = rSges(5,2) * t58 + (-t25 - t16) * t27 + t68 * t29;
t38 = cos(qJ(1));
t31 = t38 * pkin(1);
t1 = [-m(2) * (g(1) * (-t36 * rSges(2,1) - t38 * rSges(2,2)) + g(2) * (t38 * rSges(2,1) - t36 * rSges(2,2))) - m(3) * (g(1) * (t50 - t59) + g(2) * (t31 + t53)) - m(4) * (g(1) * (t45 - t59) + g(2) * (t31 + t47)) - m(5) * (g(1) * (t42 - t59) + g(2) * (t31 + t44)) - m(6) * (g(1) * (t46 - t59) + g(2) * (t31 + t43)), -m(3) * (g(1) * t50 + g(2) * t53) - m(4) * (g(1) * t45 + g(2) * t47) - m(5) * (g(1) * t42 + g(2) * t44) - m(6) * (g(1) * t46 + g(2) * t43), (-m(4) * t67 - m(5) * (t30 + t52) - m(6) * (t30 + t51)) * g(3) + t64 * (-m(4) * (-rSges(4,1) * t35 - rSges(4,2) * t37) - m(5) * (t49 - t60) - m(6) * (-t60 + t65)), (-m(5) * t52 - m(6) * t51) * g(3) + t64 * (-m(5) * t49 - m(6) * t65), -m(6) * (g(1) * t27 - g(2) * t29)];
taug = t1(:);

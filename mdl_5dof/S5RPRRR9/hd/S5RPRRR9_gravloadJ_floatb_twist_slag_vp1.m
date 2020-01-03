% Calculate Gravitation load on the joints for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:18
% EndTime: 2019-12-31 19:07:20
% DurationCPUTime: 0.42s
% Computational Cost: add. (284->78), mult. (261->104), div. (0->0), fcn. (226->10), ass. (0->43)
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t64 = g(1) * t35 + g(2) * t33;
t68 = rSges(6,3) + pkin(8);
t28 = pkin(9) + qJ(3);
t25 = qJ(4) + t28;
t20 = sin(t25);
t21 = cos(t25);
t65 = t21 * rSges(5,1) - t20 * rSges(5,2);
t39 = t21 * pkin(4) + t68 * t20;
t34 = cos(qJ(5));
t59 = rSges(6,1) * t34;
t63 = (-pkin(4) - t59) * t20;
t23 = sin(t28);
t62 = pkin(3) * t23;
t30 = cos(pkin(9));
t22 = t30 * pkin(2) + pkin(1);
t32 = sin(qJ(5));
t58 = rSges(6,2) * t32;
t54 = t33 * t32;
t53 = t33 * t34;
t52 = t35 * t32;
t51 = t35 * t34;
t31 = -pkin(6) - qJ(2);
t50 = rSges(4,3) - t31;
t27 = -pkin(7) + t31;
t49 = rSges(5,3) - t27;
t48 = rSges(3,3) + qJ(2);
t24 = cos(t28);
t44 = t24 * rSges(4,1) - t23 * rSges(4,2);
t42 = -rSges(5,1) * t20 - rSges(5,2) * t21;
t41 = rSges(3,1) * t30 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t40 = t22 + t44;
t38 = t39 + (-t58 + t59) * t21;
t37 = t64 * (t20 * t58 + t21 * t68);
t19 = pkin(3) * t24;
t8 = t19 + t22;
t5 = t35 * t8;
t4 = t21 * t51 + t54;
t3 = -t21 * t52 + t53;
t2 = -t21 * t53 + t52;
t1 = t21 * t54 + t51;
t6 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - t35 * rSges(2,2)) + g(2) * (t35 * rSges(2,1) - t33 * rSges(2,2))) - m(3) * ((g(1) * t48 + g(2) * t41) * t35 + (-g(1) * t41 + g(2) * t48) * t33) - m(4) * ((g(1) * t50 + g(2) * t40) * t35 + (-g(1) * t40 + g(2) * t50) * t33) - m(5) * (g(2) * t5 + (g(1) * t49 + g(2) * t65) * t35 + (g(1) * (-t65 - t8) + g(2) * t49) * t33) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2)) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t5) + (-g(1) * t27 + g(2) * t39) * t35 + (g(1) * (-t39 - t8) - g(2) * t27) * t33), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t33 - g(2) * t35), -m(6) * t37 + (-m(4) * t44 - m(5) * (t19 + t65) - m(6) * (t19 + t38)) * g(3) + t64 * (-m(4) * (-rSges(4,1) * t23 - rSges(4,2) * t24) - m(5) * (t42 - t62) - m(6) * (-t62 + t63)), -m(5) * (g(3) * t65 + t64 * t42) - m(6) * (g(3) * t38 + t64 * t63 + t37), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t32 - rSges(6,2) * t34) * t20)];
taug = t6(:);

% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:58:55
% DurationCPUTime: 0.24s
% Computational Cost: add. (244->63), mult. (166->77), div. (0->0), fcn. (127->8), ass. (0->32)
t56 = rSges(6,1) + pkin(4);
t34 = cos(qJ(4));
t55 = t56 * t34;
t54 = rSges(5,3) + pkin(7);
t49 = rSges(5,1) * t34;
t53 = pkin(3) + t49;
t52 = rSges(6,3) + qJ(5) + pkin(7);
t51 = pkin(3) + t55;
t30 = qJ(1) + qJ(2);
t24 = pkin(8) + t30;
t19 = sin(t24);
t32 = sin(qJ(4));
t48 = t19 * t32;
t20 = cos(t24);
t47 = t20 * t32;
t25 = sin(t30);
t26 = cos(t30);
t45 = t25 * rSges(3,1) + t26 * rSges(3,2);
t21 = pkin(2) * t25;
t44 = t19 * rSges(4,1) + t20 * rSges(4,2) + t21;
t43 = t26 * rSges(3,1) - rSges(3,2) * t25;
t22 = pkin(2) * t26;
t42 = t20 * rSges(4,1) - rSges(4,2) * t19 + t22;
t39 = -rSges(5,2) * t47 + t54 * t19 + t53 * t20 + t22;
t38 = -rSges(6,2) * t48 + t51 * t19 - t52 * t20 + t21;
t37 = -rSges(6,2) * t47 + t52 * t19 + t51 * t20 + t22;
t36 = -rSges(5,2) * t48 + t53 * t19 - t54 * t20 + t21;
t35 = cos(qJ(1));
t33 = sin(qJ(1));
t29 = t35 * pkin(1);
t27 = t33 * pkin(1);
t1 = [-m(2) * (g(2) * (rSges(2,1) * t35 - t33 * rSges(2,2)) + g(3) * (t33 * rSges(2,1) + rSges(2,2) * t35)) - m(3) * (g(2) * (t29 + t43) + g(3) * (t27 + t45)) - m(4) * (g(2) * (t29 + t42) + g(3) * (t27 + t44)) - m(5) * (g(2) * (t29 + t39) + g(3) * (t27 + t36)) - m(6) * (g(2) * (t29 + t37) + g(3) * (t27 + t38)), -m(3) * (g(2) * t43 + g(3) * t45) - m(4) * (g(2) * t42 + g(3) * t44) - m(5) * (g(2) * t39 + g(3) * t36) - m(6) * (g(2) * t37 + g(3) * t38), (-m(4) - m(5) - m(6)) * g(1), (-m(5) * (-rSges(5,2) * t32 + t49) - m(6) * (-rSges(6,2) * t32 + t55)) * g(1) + (g(2) * t19 - g(3) * t20) * (m(5) * (rSges(5,1) * t32 + rSges(5,2) * t34) + m(6) * (rSges(6,2) * t34 + t56 * t32)), -m(6) * (-g(2) * t20 - g(3) * t19)];
taug = t1(:);

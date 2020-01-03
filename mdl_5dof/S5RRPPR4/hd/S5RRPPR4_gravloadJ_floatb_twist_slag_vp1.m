% Calculate Gravitation load on the joints for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:41
% DurationCPUTime: 0.18s
% Computational Cost: add. (248->58), mult. (212->71), div. (0->0), fcn. (208->8), ass. (0->30)
t53 = rSges(6,3) + pkin(7);
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t52 = -rSges(6,1) * t31 + rSges(6,2) * t29;
t51 = pkin(4) - t52;
t50 = m(5) + m(6);
t30 = sin(qJ(1));
t49 = pkin(1) * t30;
t28 = qJ(1) + qJ(2);
t25 = sin(t28);
t26 = cos(t28);
t46 = t26 * pkin(2) + t25 * qJ(3);
t45 = cos(pkin(8));
t44 = sin(pkin(8));
t43 = t26 * pkin(3) + t46;
t42 = t26 * rSges(3,1) - rSges(3,2) * t25;
t41 = t26 * rSges(4,1) + t25 * rSges(4,3) + t46;
t11 = -t25 * t44 - t26 * t45;
t12 = -t25 * t45 + t26 * t44;
t40 = -t11 * rSges(5,1) - t12 * rSges(5,2) + t43;
t18 = t26 * qJ(3);
t39 = t18 + (-pkin(2) - pkin(3)) * t25;
t38 = -rSges(3,1) * t25 - rSges(3,2) * t26;
t36 = t18 + t26 * rSges(4,3) + (-rSges(4,1) - pkin(2)) * t25;
t35 = t12 * rSges(5,1) - t11 * rSges(5,2) + t39;
t34 = -t51 * t11 + t53 * t12 + t43;
t33 = t53 * t11 + t51 * t12 + t39;
t32 = cos(qJ(1));
t27 = t32 * pkin(1);
t1 = [-m(2) * (g(1) * (-t30 * rSges(2,1) - rSges(2,2) * t32) + g(2) * (rSges(2,1) * t32 - t30 * rSges(2,2))) - m(3) * (g(1) * (t38 - t49) + g(2) * (t27 + t42)) - m(4) * (g(1) * (t36 - t49) + g(2) * (t27 + t41)) - m(5) * (g(1) * (t35 - t49) + g(2) * (t27 + t40)) - m(6) * (g(1) * (t33 - t49) + g(2) * (t27 + t34)), -m(3) * (g(1) * t38 + g(2) * t42) - m(4) * (g(1) * t36 + g(2) * t41) - m(5) * (g(1) * t35 + g(2) * t40) - m(6) * (g(1) * t33 + g(2) * t34), (-m(4) - t50) * (g(1) * t25 - g(2) * t26), t50 * g(3), -m(6) * (g(3) * t52 + (g(1) * t11 + g(2) * t12) * (rSges(6,1) * t29 + rSges(6,2) * t31))];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (215->68), mult. (196->81), div. (0->0), fcn. (157->6), ass. (0->33)
t54 = rSges(6,1) + pkin(4);
t23 = qJ(1) + qJ(2);
t21 = cos(t23);
t20 = sin(t23);
t49 = g(1) * t20;
t55 = g(2) * t21 - t49;
t26 = cos(qJ(4));
t42 = rSges(6,3) + qJ(5);
t53 = t42 * t26;
t51 = -pkin(2) - pkin(7);
t25 = sin(qJ(1));
t50 = pkin(1) * t25;
t47 = rSges(5,2) * t26;
t24 = sin(qJ(4));
t46 = t20 * t24;
t45 = t21 * t24;
t44 = t21 * t26;
t43 = t21 * pkin(2) + t20 * qJ(3);
t11 = t21 * qJ(3);
t41 = rSges(5,1) * t45 + rSges(5,2) * t44 + t11;
t40 = t21 * pkin(7) + t43;
t39 = t21 * rSges(3,1) - rSges(3,2) * t20;
t37 = (-rSges(5,3) + t51) * t49;
t36 = -rSges(3,1) * t20 - rSges(3,2) * t21;
t34 = t21 * rSges(4,3) + t11 + (rSges(4,2) - pkin(2)) * t20;
t33 = t21 * rSges(6,2) + t54 * t46 + t40;
t32 = rSges(5,1) * t46 + t21 * rSges(5,3) + t20 * t47 + t40;
t31 = -rSges(4,2) * t21 + t20 * rSges(4,3) + t43;
t30 = -t42 * t44 + t54 * t45 + t11;
t28 = (g(1) * (-rSges(6,2) + t51) - g(2) * t53) * t20;
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t25 * rSges(2,2))) - m(3) * (g(1) * (t36 - t50) + g(2) * (t22 + t39)) - m(4) * (g(1) * (t34 - t50) + g(2) * (t22 + t31)) - m(5) * (g(1) * (t41 - t50) + g(2) * (t22 + t32) + t37) - m(6) * (g(1) * (t30 - t50) + g(2) * (t22 + t33) + t28), -m(3) * (g(1) * t36 + g(2) * t39) - m(4) * (g(1) * t34 + g(2) * t31) - m(5) * (g(1) * t41 + g(2) * t32 + t37) - m(6) * (g(1) * t30 + g(2) * t33 + t28), -(-m(4) - m(5) - m(6)) * t55, (-m(5) * (-rSges(5,1) * t24 - t47) - m(6) * (-t54 * t24 + t53)) * g(3) + t55 * (m(5) * (rSges(5,1) * t26 - rSges(5,2) * t24) + m(6) * (t42 * t24 + t54 * t26)), -m(6) * (g(3) * t24 + t26 * t55)];
taug = t1(:);

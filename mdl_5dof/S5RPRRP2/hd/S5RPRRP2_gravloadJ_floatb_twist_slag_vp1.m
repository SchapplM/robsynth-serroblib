% Calculate Gravitation load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:07
% DurationCPUTime: 0.20s
% Computational Cost: add. (225->62), mult. (154->81), div. (0->0), fcn. (117->8), ass. (0->29)
t49 = rSges(5,3) + pkin(7);
t48 = rSges(6,3) + qJ(5) + pkin(7);
t27 = qJ(1) + pkin(8);
t23 = qJ(3) + t27;
t18 = sin(t23);
t29 = sin(qJ(4));
t46 = t18 * t29;
t31 = cos(qJ(4));
t45 = t18 * t31;
t19 = cos(t23);
t44 = t19 * t29;
t43 = t19 * t31;
t42 = t18 * rSges(4,1) + t19 * rSges(4,2);
t21 = sin(t27);
t30 = sin(qJ(1));
t24 = t30 * pkin(1);
t41 = pkin(2) * t21 + t24;
t22 = cos(t27);
t32 = cos(qJ(1));
t26 = t32 * pkin(1);
t40 = pkin(2) * t22 + t26;
t39 = t19 * rSges(4,1) - rSges(4,2) * t18;
t37 = rSges(5,1) * t43 - rSges(5,2) * t44 + t19 * pkin(3) + t49 * t18;
t25 = t31 * pkin(4);
t20 = t25 + pkin(3);
t35 = rSges(6,1) * t45 - rSges(6,2) * t46 + t18 * t20 - t48 * t19;
t34 = rSges(6,1) * t43 - rSges(6,2) * t44 + t48 * t18 + t19 * t20;
t33 = rSges(5,1) * t45 - rSges(5,2) * t46 + t18 * pkin(3) - t49 * t19;
t1 = [-m(2) * (g(2) * (rSges(2,1) * t32 - t30 * rSges(2,2)) + g(3) * (t30 * rSges(2,1) + rSges(2,2) * t32)) - m(3) * (g(2) * (rSges(3,1) * t22 - rSges(3,2) * t21 + t26) + g(3) * (rSges(3,1) * t21 + rSges(3,2) * t22 + t24)) - m(4) * (g(2) * (t39 + t40) + g(3) * (t41 + t42)) - m(5) * (g(2) * (t37 + t40) + g(3) * (t33 + t41)) - m(6) * (g(2) * (t34 + t40) + g(3) * (t35 + t41)), (-m(3) - m(4) - m(5) - m(6)) * g(1), -m(4) * (g(2) * t39 + g(3) * t42) - m(5) * (g(2) * t37 + g(3) * t33) - m(6) * (g(2) * t34 + g(3) * t35), (-m(5) * (rSges(5,1) * t31 - rSges(5,2) * t29) - m(6) * (rSges(6,1) * t31 - rSges(6,2) * t29 + t25)) * g(1) + (g(2) * t18 - g(3) * t19) * (m(5) * (rSges(5,1) * t29 + rSges(5,2) * t31) + m(6) * (rSges(6,2) * t31 + (rSges(6,1) + pkin(4)) * t29)), -m(6) * (-g(2) * t19 - g(3) * t18)];
taug = t1(:);

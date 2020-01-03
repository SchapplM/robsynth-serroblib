% Calculate Gravitation load on the joints for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:38
% EndTime: 2019-12-31 17:35:39
% DurationCPUTime: 0.15s
% Computational Cost: add. (155->33), mult. (174->50), div. (0->0), fcn. (182->8), ass. (0->24)
t44 = -rSges(6,3) - pkin(7);
t20 = cos(pkin(8));
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t35 = sin(pkin(8));
t43 = -t20 * t22 + t35 * t24;
t21 = sin(qJ(5));
t23 = cos(qJ(5));
t42 = -rSges(6,1) * t23 + rSges(6,2) * t21;
t41 = pkin(4) - t42;
t34 = qJ(3) + qJ(4);
t31 = sin(t34);
t32 = cos(t34);
t11 = -t20 * t32 - t35 * t31;
t12 = t20 * t31 - t35 * t32;
t40 = t11 * rSges(5,1) + t12 * rSges(5,2);
t39 = -t12 * rSges(5,1) + t11 * rSges(5,2);
t33 = -m(3) - m(4) - m(5) - m(6);
t29 = t43 * pkin(3);
t13 = -t20 * t24 - t35 * t22;
t27 = t13 * pkin(3);
t26 = t41 * t11 + t44 * t12;
t25 = t44 * t11 - t41 * t12;
t1 = [(-m(2) + t33) * g(3), t33 * (g(1) * t35 - g(2) * t20), -m(4) * (g(1) * (rSges(4,1) * t43 + rSges(4,2) * t13) + g(2) * (rSges(4,1) * t13 - rSges(4,2) * t43)) - m(5) * (g(1) * (t29 + t39) + g(2) * (t27 + t40)) - m(6) * (g(1) * (t25 + t29) + g(2) * (t27 + t26)), -m(5) * (g(1) * t39 + g(2) * t40) - m(6) * (g(1) * t25 + g(2) * t26), -m(6) * (g(3) * t42 + (g(1) * t11 + g(2) * t12) * (rSges(6,1) * t21 + rSges(6,2) * t23))];
taug = t1(:);

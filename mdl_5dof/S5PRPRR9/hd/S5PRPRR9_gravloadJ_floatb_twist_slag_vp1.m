% Calculate Gravitation load on the joints for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (191->41), mult. (170->57), div. (0->0), fcn. (172->6), ass. (0->21)
t23 = pkin(8) + qJ(2);
t11 = cos(t23);
t21 = sin(t23);
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t1 = -t11 * t26 - t21 * t25;
t12 = sin(qJ(5));
t13 = cos(qJ(5));
t18 = -rSges(6,1) * t13 + t12 * rSges(6,2);
t16 = pkin(4) - t18;
t2 = t11 * t25 - t21 * t26;
t28 = rSges(6,3) + pkin(7);
t31 = -(g(1) * t28 - g(2) * t16) * t1 - (g(1) * t16 + g(2) * t28) * t2;
t27 = t11 * pkin(2) + t21 * qJ(3);
t24 = -m(4) - m(5) - m(6);
t22 = t11 * pkin(3) + t27;
t20 = -t2 * rSges(5,1) + t1 * rSges(5,2);
t19 = rSges(5,1) * t1 + rSges(5,2) * t2;
t15 = -t21 * pkin(2) + t11 * qJ(3);
t14 = -t21 * pkin(3) + t15;
t3 = [(-m(2) - m(3) + t24) * g(3), -m(3) * (g(1) * (-t21 * rSges(3,1) - t11 * rSges(3,2)) + g(2) * (t11 * rSges(3,1) - t21 * rSges(3,2))) - m(4) * (g(1) * (-t21 * rSges(4,1) + t11 * rSges(4,3) + t15) + g(2) * (t11 * rSges(4,1) + t21 * rSges(4,3) + t27)) - m(5) * (g(1) * (t14 - t20) + g(2) * (-t19 + t22)) - m(6) * (g(1) * t14 + g(2) * t22 - t31), t24 * (g(1) * t21 - g(2) * t11), -m(5) * (g(1) * t20 + g(2) * t19) - m(6) * t31, -m(6) * (g(3) * t18 + (g(1) * t1 + g(2) * t2) * (rSges(6,1) * t12 + rSges(6,2) * t13))];
taug = t3(:);

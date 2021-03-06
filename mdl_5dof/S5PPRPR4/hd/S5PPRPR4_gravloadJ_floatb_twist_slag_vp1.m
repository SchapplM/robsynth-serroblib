% Calculate Gravitation load on the joints for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:16
% EndTime: 2019-12-31 17:32:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (89->30), mult. (144->44), div. (0->0), fcn. (152->8), ass. (0->18)
t26 = -m(5) - m(6);
t25 = cos(qJ(3));
t24 = sin(qJ(3));
t23 = -rSges(6,3) - pkin(6) - qJ(4);
t22 = -rSges(5,3) - qJ(4);
t21 = cos(pkin(7));
t20 = sin(pkin(7));
t19 = -m(3) - m(4) + t26;
t11 = pkin(8) + qJ(5);
t10 = cos(t11);
t9 = sin(t11);
t17 = -rSges(6,1) * t10 + rSges(6,2) * t9;
t13 = cos(pkin(8));
t16 = pkin(4) * t13 + pkin(3) - t17;
t15 = rSges(5,1) * t13 - rSges(5,2) * sin(pkin(8)) + pkin(3);
t3 = -t20 * t25 + t21 * t24;
t2 = -t20 * t24 - t21 * t25;
t1 = [(-m(2) + t19) * g(3), t19 * (g(1) * t20 - g(2) * t21), -m(4) * (g(1) * (-rSges(4,1) * t3 + rSges(4,2) * t2) + g(2) * (rSges(4,1) * t2 + rSges(4,2) * t3)) - m(5) * ((-g(1) * t15 + g(2) * t22) * t3 + (g(1) * t22 + g(2) * t15) * t2) - m(6) * ((-g(1) * t16 + g(2) * t23) * t3 + (g(1) * t23 + g(2) * t16) * t2), t26 * (g(1) * t3 - g(2) * t2), -m(6) * (g(3) * t17 + (g(1) * t2 + g(2) * t3) * (rSges(6,1) * t9 + rSges(6,2) * t10))];
taug = t1(:);

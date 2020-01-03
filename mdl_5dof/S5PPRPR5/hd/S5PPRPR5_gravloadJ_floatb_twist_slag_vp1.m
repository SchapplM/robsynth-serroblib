% Calculate Gravitation load on the joints for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:24
% EndTime: 2019-12-31 17:33:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (75->29), mult. (138->43), div. (0->0), fcn. (146->6), ass. (0->18)
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t19 = sin(qJ(3));
t20 = cos(qJ(3));
t4 = -t16 * t19 - t17 * t20;
t5 = -t16 * t20 + t17 * t19;
t25 = g(1) * t5 - g(2) * t4;
t22 = -m(5) - m(6);
t21 = rSges(6,3) + pkin(6);
t18 = -rSges(5,3) - qJ(4);
t15 = -m(3) - m(4) + t22;
t10 = sin(qJ(5));
t11 = cos(qJ(5));
t13 = t10 * rSges(6,1) + rSges(6,2) * t11;
t12 = -qJ(4) - t13;
t3 = t5 * pkin(3);
t2 = t4 * pkin(3);
t1 = [(-m(2) + t15) * g(3), t15 * (g(1) * t16 - g(2) * t17), -m(4) * (g(1) * (-rSges(4,1) * t5 + rSges(4,2) * t4) + g(2) * (rSges(4,1) * t4 + rSges(4,2) * t5)) - m(5) * (g(1) * (rSges(5,2) * t5 + t18 * t4 - t3) + g(2) * (-rSges(5,2) * t4 + t18 * t5 + t2)) - m(6) * (-g(1) * t3 + g(2) * t2 + (-g(1) * t21 + g(2) * t12) * t5 + (g(1) * t12 + g(2) * t21) * t4), t22 * t25, -m(6) * (g(3) * t13 + t25 * (rSges(6,1) * t11 - rSges(6,2) * t10))];
taug = t1(:);

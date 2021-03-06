% Calculate Gravitation load on the joints for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (90->31), mult. (168->48), div. (0->0), fcn. (177->6), ass. (0->17)
t10 = cos(qJ(4));
t23 = rSges(6,1) + pkin(4);
t9 = sin(qJ(4));
t24 = t9 * rSges(6,2) - t23 * t10;
t21 = -rSges(5,3) - pkin(6);
t20 = -rSges(6,3) - qJ(5) - pkin(6);
t19 = cos(qJ(3));
t18 = sin(qJ(3));
t17 = cos(pkin(7));
t16 = sin(pkin(7));
t15 = -m(3) - m(4) - m(5) - m(6);
t14 = -rSges(5,1) * t10 + t9 * rSges(5,2);
t13 = pkin(3) - t14;
t12 = pkin(3) - t24;
t2 = -t16 * t19 + t17 * t18;
t1 = -t16 * t18 - t17 * t19;
t3 = [(-m(2) + t15) * g(3), t15 * (g(1) * t16 - g(2) * t17), -m(4) * (g(1) * (-rSges(4,1) * t2 + rSges(4,2) * t1) + g(2) * (rSges(4,1) * t1 + rSges(4,2) * t2)) - m(5) * ((-g(1) * t13 + g(2) * t21) * t2 + (g(1) * t21 + g(2) * t13) * t1) - m(6) * ((-g(1) * t12 + g(2) * t20) * t2 + (g(1) * t20 + g(2) * t12) * t1), (-m(5) * t14 - m(6) * t24) * g(3) + (g(1) * t1 + g(2) * t2) * (-m(5) * (rSges(5,1) * t9 + rSges(5,2) * t10) - m(6) * (rSges(6,2) * t10 + t23 * t9)), -m(6) * (g(1) * t2 - g(2) * t1)];
taug = t3(:);

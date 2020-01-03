% Calculate Gravitation load on the joints for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:03
% DurationCPUTime: 0.22s
% Computational Cost: add. (64->42), mult. (117->55), div. (0->0), fcn. (93->4), ass. (0->14)
t7 = sin(qJ(1));
t9 = cos(qJ(1));
t26 = -g(1) * t7 + g(2) * t9;
t23 = rSges(5,3) + qJ(4);
t24 = rSges(5,1) + pkin(3);
t6 = sin(qJ(3));
t8 = cos(qJ(3));
t25 = t23 * t8 - t24 * t6;
t12 = rSges(4,1) * t6 + rSges(4,2) * t8;
t19 = -pkin(1) - pkin(5);
t18 = t9 * pkin(1) + t7 * qJ(2);
t15 = t9 * pkin(5) + t18;
t3 = t9 * qJ(2);
t1 = [-m(2) * (g(1) * (-t7 * rSges(2,1) - rSges(2,2) * t9) + g(2) * (rSges(2,1) * t9 - t7 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t9 + t3 + (rSges(3,2) - pkin(1)) * t7) + g(2) * (-rSges(3,2) * t9 + t7 * rSges(3,3) + t18)) - m(4) * (g(1) * (t12 * t9 + t3) + g(2) * (rSges(4,3) * t9 + t15) + (g(1) * (-rSges(4,3) + t19) + g(2) * t12) * t7) - m(5) * (g(1) * t3 + g(2) * t15 + (g(2) * rSges(5,2) - g(1) * t25) * t9 + (g(1) * (-rSges(5,2) + t19) - g(2) * t25) * t7), -(-m(3) - m(4) - m(5)) * t26, (m(4) * t12 - m(5) * t25) * g(3) + t26 * (m(4) * (rSges(4,1) * t8 - rSges(4,2) * t6) + m(5) * (t23 * t6 + t24 * t8)), -m(5) * (g(3) * t6 + t26 * t8)];
taug = t1(:);

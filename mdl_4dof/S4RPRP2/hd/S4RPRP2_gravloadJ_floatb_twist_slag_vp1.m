% Calculate Gravitation load on the joints for
% S4RPRP2
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:42
% EndTime: 2019-03-08 18:30:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (67->42), mult. (114->55), div. (0->0), fcn. (108->4), ass. (0->16)
t12 = sin(qJ(1));
t13 = cos(qJ(1));
t21 = t13 * pkin(1) + t12 * qJ(2);
t20 = cos(qJ(3));
t11 = sin(qJ(3));
t19 = t12 * t11;
t18 = t13 * t11;
t1 = -t13 * t20 - t19;
t2 = -t12 * t20 + t18;
t17 = -rSges(4,1) * t2 + rSges(4,2) * t1;
t16 = t1 * rSges(4,1) + t2 * rSges(4,2);
t15 = -t2 * rSges(5,1) + t1 * rSges(5,2);
t14 = t1 * rSges(5,1) + t2 * rSges(5,2);
t9 = t13 * qJ(2);
t7 = t20 * pkin(3) + pkin(2);
t3 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - rSges(2,2) * t13) + g(2) * (rSges(2,1) * t13 - t12 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t13 + t9 + (-rSges(3,1) - pkin(1)) * t12) + g(2) * (rSges(3,1) * t13 + t12 * rSges(3,3) + t21)) - m(4) * (g(1) * (t9 + (-pkin(1) - pkin(2)) * t12 - t17) + g(2) * (pkin(2) * t13 - t16 + t21)) - m(5) * (g(1) * (pkin(3) * t18 + t9 + (-pkin(1) - t7) * t12 - t15) + g(2) * (pkin(3) * t19 + t13 * t7 - t14 + t21)) (-m(3) - m(4) - m(5)) * (g(1) * t12 - g(2) * t13) -m(4) * (g(1) * t17 + g(2) * t16) - m(5) * (g(1) * t15 + g(2) * t14 + (-g(1) * t2 + g(2) * t1) * pkin(3)) m(5) * g(3)];
taug  = t3(:);

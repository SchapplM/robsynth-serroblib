% Calculate Gravitation load on the joints for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:32
% EndTime: 2019-12-31 16:51:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (93->39), mult. (165->56), div. (0->0), fcn. (172->6), ass. (0->19)
t13 = cos(qJ(1));
t22 = sin(qJ(3));
t23 = sin(qJ(1));
t24 = cos(qJ(3));
t1 = -t13 * t24 - t23 * t22;
t11 = sin(qJ(4));
t12 = cos(qJ(4));
t17 = -rSges(5,1) * t12 + rSges(5,2) * t11;
t15 = pkin(3) - t17;
t2 = t13 * t22 - t23 * t24;
t26 = rSges(5,3) + pkin(6);
t29 = -(g(1) * t26 - g(2) * t15) * t1 - (g(1) * t15 + g(2) * t26) * t2;
t25 = t13 * pkin(1) + t23 * qJ(2);
t21 = t13 * pkin(2) + t25;
t20 = -t23 * pkin(1) + t13 * qJ(2);
t19 = -t2 * rSges(4,1) + t1 * rSges(4,2);
t18 = rSges(4,1) * t1 + rSges(4,2) * t2;
t14 = -t23 * pkin(2) + t20;
t3 = [-m(2) * (g(1) * (-t23 * rSges(2,1) - t13 * rSges(2,2)) + g(2) * (t13 * rSges(2,1) - t23 * rSges(2,2))) - m(3) * (g(1) * (-t23 * rSges(3,1) + t13 * rSges(3,3) + t20) + g(2) * (t13 * rSges(3,1) + t23 * rSges(3,3) + t25)) - m(4) * (g(1) * (t14 - t19) + g(2) * (-t18 + t21)) - m(5) * (g(1) * t14 + g(2) * t21 - t29), (-m(3) - m(4) - m(5)) * (g(1) * t23 - g(2) * t13), -m(4) * (g(1) * t19 + g(2) * t18) - m(5) * t29, -m(5) * (g(3) * t17 + (g(1) * t1 + g(2) * t2) * (rSges(5,1) * t11 + rSges(5,2) * t12))];
taug = t3(:);

% Calculate Gravitation load on the joints for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:38
% EndTime: 2019-03-08 18:15:38
% DurationCPUTime: 0.09s
% Computational Cost: add. (51->18), mult. (66->30), div. (0->0), fcn. (62->6), ass. (0->14)
t19 = -m(3) - m(4) - m(5);
t13 = sin(pkin(6));
t14 = cos(pkin(6));
t15 = sin(qJ(3));
t16 = cos(qJ(3));
t18 = t13 * t16 - t14 * t15;
t7 = -t13 * t15 - t14 * t16;
t12 = qJ(3) + qJ(4);
t10 = sin(t12);
t11 = cos(t12);
t5 = -t13 * t10 - t14 * t11;
t6 = t14 * t10 - t13 * t11;
t17 = g(1) * (-t6 * rSges(5,1) + t5 * rSges(5,2)) + g(2) * (t5 * rSges(5,1) + t6 * rSges(5,2));
t1 = [(-m(2) + t19) * g(3), t19 * (g(1) * t13 - g(2) * t14) -m(4) * (g(1) * (rSges(4,1) * t18 + t7 * rSges(4,2)) + g(2) * (t7 * rSges(4,1) - rSges(4,2) * t18)) - m(5) * ((g(1) * t18 + g(2) * t7) * pkin(3) + t17) -m(5) * t17];
taug  = t1(:);

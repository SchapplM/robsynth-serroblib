% Calculate Gravitation load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:20
% EndTime: 2019-03-08 18:26:20
% DurationCPUTime: 0.17s
% Computational Cost: add. (89->48), mult. (189->68), div. (0->0), fcn. (199->6), ass. (0->26)
t35 = rSges(5,2) + qJ(3);
t34 = rSges(5,3) + qJ(4);
t33 = -m(4) - m(5);
t15 = sin(pkin(4));
t32 = g(3) * t15;
t14 = sin(pkin(6));
t18 = sin(qJ(1));
t30 = t14 * t18;
t29 = t15 * t18;
t19 = cos(qJ(1));
t28 = t15 * t19;
t17 = cos(pkin(4));
t27 = t17 * t19;
t16 = cos(pkin(6));
t26 = t18 * t16;
t24 = qJ(2) * t15;
t25 = t19 * pkin(1) + t18 * t24;
t23 = rSges(4,3) + qJ(3);
t8 = t16 * t19 - t17 * t30;
t22 = t8 * pkin(2) + t25;
t21 = -t18 * pkin(1) + t19 * t24;
t6 = t14 * t27 + t26;
t20 = -t6 * pkin(2) + t21;
t7 = t14 * t19 + t17 * t26;
t5 = -t16 * t27 + t30;
t1 = [-m(2) * (g(1) * (-t18 * rSges(2,1) - rSges(2,2) * t19) + g(2) * (rSges(2,1) * t19 - t18 * rSges(2,2))) - m(3) * (g(1) * (-t6 * rSges(3,1) + t5 * rSges(3,2) + rSges(3,3) * t28 + t21) + g(2) * (rSges(3,1) * t8 - rSges(3,2) * t7 + rSges(3,3) * t29 + t25)) - m(4) * (g(1) * (rSges(4,1) * t28 + t6 * rSges(4,2) - t23 * t5 + t20) + g(2) * (rSges(4,1) * t29 - rSges(4,2) * t8 + t23 * t7 + t22)) - m(5) * (g(1) * (-t34 * t6 - t35 * t5 + t20) + g(2) * (t34 * t8 + t35 * t7 + t22) + (g(1) * t19 + g(2) * t18) * t15 * (rSges(5,1) + pkin(3))) (-m(3) + t33) * (g(3) * t17 + (g(1) * t18 - g(2) * t19) * t15) t33 * (g(1) * t7 + g(2) * t5 - t16 * t32) -m(5) * (g(1) * t8 + g(2) * t6 + t14 * t32)];
taug  = t1(:);

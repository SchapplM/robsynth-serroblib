% Calculate Gravitation load on the joints for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:37
% DurationCPUTime: 0.17s
% Computational Cost: add. (53->35), mult. (118->52), div. (0->0), fcn. (100->6), ass. (0->16)
t6 = sin(pkin(6));
t7 = cos(pkin(6));
t24 = g(1) * t7 + g(2) * t6;
t8 = sin(qJ(4));
t9 = sin(qJ(2));
t21 = t8 * t9;
t20 = -m(4) - m(5);
t11 = cos(qJ(2));
t19 = t11 * pkin(2) + t9 * qJ(3);
t18 = g(3) * t11;
t10 = cos(qJ(4));
t17 = t10 * t9;
t16 = rSges(5,3) + pkin(5);
t14 = t24 * qJ(3) * t11;
t13 = rSges(5,1) * t8 + rSges(5,2) * t10;
t1 = [(-m(2) - m(3) + t20) * g(3), -m(3) * g(3) * (rSges(3,1) * t11 - t9 * rSges(3,2)) - m(4) * (g(3) * (-rSges(4,2) * t11 + t9 * rSges(4,3) + t19) + t14) - m(5) * (g(3) * (t16 * t11 + t13 * t9 + t19) + t14) + t24 * ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * t13) * t11 + (m(3) * rSges(3,1) - m(4) * (rSges(4,2) - pkin(2)) - m(5) * (-pkin(2) - t16)) * t9), t20 * (t24 * t9 - t18), -m(5) * (g(1) * ((t7 * t17 - t6 * t8) * rSges(5,1) + (-t10 * t6 - t7 * t21) * rSges(5,2)) + g(2) * ((t6 * t17 + t7 * t8) * rSges(5,1) + (t10 * t7 - t6 * t21) * rSges(5,2)) + (-rSges(5,1) * t10 + rSges(5,2) * t8) * t18)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:58:54
% DurationCPUTime: 0.28s
% Computational Cost: add. (84->57), mult. (162->74), div. (0->0), fcn. (135->4), ass. (0->23)
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t32 = g(1) * t14 + g(2) * t12;
t13 = cos(qJ(2));
t11 = sin(qJ(2));
t26 = rSges(5,2) * t11;
t28 = rSges(5,1) + pkin(3);
t33 = t28 * t13 + t26;
t6 = t11 * qJ(3);
t8 = t13 * pkin(2);
t31 = t6 + t8;
t27 = t14 * pkin(1) + t12 * pkin(5);
t25 = t11 * t14;
t24 = t13 * t14;
t22 = -rSges(5,3) - qJ(4);
t21 = -pkin(2) - t28;
t20 = pkin(2) * t24 + t14 * t6 + t27;
t19 = -pkin(1) - t6;
t18 = t32 * qJ(3) * t13;
t17 = rSges(3,1) * t13 - rSges(3,2) * t11;
t16 = rSges(4,1) * t13 + rSges(4,3) * t11;
t9 = t14 * pkin(5);
t1 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - rSges(2,2) * t14) + g(2) * (rSges(2,1) * t14 - t12 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t14 + t9) + g(2) * (rSges(3,1) * t24 - rSges(3,2) * t25 + t27) + (g(1) * (-pkin(1) - t17) + g(2) * rSges(3,3)) * t12) - m(4) * (g(1) * (rSges(4,2) * t14 + t9) + g(2) * (rSges(4,1) * t24 + rSges(4,3) * t25 + t20) + (g(1) * (-t16 + t19 - t8) + g(2) * rSges(4,2)) * t12) - m(5) * (g(1) * t9 + g(2) * t20 + (g(1) * t22 + g(2) * t33) * t14 + (g(2) * t22 + (t21 * t13 + t19 - t26) * g(1)) * t12), -m(3) * g(3) * t17 - m(4) * (g(3) * (t16 + t31) + t18) - m(5) * (g(3) * (t31 + t33) + t18) + t32 * ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * rSges(5,2)) * t13 + (m(3) * rSges(3,1) - m(4) * (-rSges(4,1) - pkin(2)) - m(5) * t21) * t11), (-m(4) - m(5)) * (-g(3) * t13 + t32 * t11), -m(5) * (-g(1) * t12 + g(2) * t14)];
taug = t1(:);

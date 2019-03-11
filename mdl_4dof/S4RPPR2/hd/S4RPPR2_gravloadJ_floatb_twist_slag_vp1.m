% Calculate Gravitation load on the joints for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:26
% EndTime: 2019-03-08 18:28:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (73->42), mult. (92->53), div. (0->0), fcn. (84->6), ass. (0->20)
t25 = m(4) + m(5);
t13 = sin(pkin(6));
t16 = cos(qJ(1));
t24 = t13 * t16;
t15 = sin(qJ(1));
t23 = t15 * t13;
t22 = t16 * pkin(1) + t15 * qJ(2);
t21 = pkin(6) + qJ(4);
t20 = cos(t21);
t19 = sin(t21);
t1 = -t15 * t19 - t16 * t20;
t2 = -t15 * t20 + t16 * t19;
t18 = -t2 * rSges(5,1) + t1 * rSges(5,2);
t17 = t1 * rSges(5,1) + t2 * rSges(5,2);
t14 = cos(pkin(6));
t11 = t16 * qJ(2);
t9 = pkin(3) * t14 + pkin(2);
t4 = t14 * t16 + t23;
t3 = -t15 * t14 + t24;
t5 = [-m(2) * (g(1) * (-t15 * rSges(2,1) - t16 * rSges(2,2)) + g(2) * (t16 * rSges(2,1) - t15 * rSges(2,2))) - m(3) * (g(1) * (t16 * rSges(3,3) + t11 + (-rSges(3,1) - pkin(1)) * t15) + g(2) * (t16 * rSges(3,1) + t15 * rSges(3,3) + t22)) - m(4) * (g(1) * (t3 * rSges(4,1) + t4 * rSges(4,2) + t11 + (-pkin(1) - pkin(2)) * t15) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t16 * pkin(2) + t22)) - m(5) * (g(1) * (pkin(3) * t24 + t11 + (-pkin(1) - t9) * t15 - t18) + g(2) * (pkin(3) * t23 + t16 * t9 - t17 + t22)) (-m(3) - t25) * (g(1) * t15 - g(2) * t16) t25 * g(3), -m(5) * (g(1) * t18 + g(2) * t17)];
taug  = t5(:);

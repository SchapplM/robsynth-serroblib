% Calculate Gravitation load on the joints for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:09
% DurationCPUTime: 0.39s
% Computational Cost: add. (181->58), mult. (330->78), div. (0->0), fcn. (363->6), ass. (0->28)
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t41 = rSges(6,3) + qJ(5);
t43 = rSges(6,1) + pkin(4);
t19 = t41 * t15 + t43 * t16;
t34 = rSges(6,2) + pkin(7);
t29 = sin(qJ(3));
t30 = sin(qJ(1));
t31 = cos(qJ(3));
t32 = cos(qJ(1));
t5 = -t30 * t29 - t32 * t31;
t6 = t32 * t29 - t30 * t31;
t49 = -(-g(1) * t34 + g(2) * t19) * t5 + (g(1) * t19 + g(2) * t34) * t6;
t22 = t16 * rSges(5,1) - t15 * rSges(5,2);
t33 = rSges(5,3) + pkin(7);
t46 = (g(1) * t22 + g(2) * t33) * t6 + (g(1) * t33 - g(2) * t22) * t5;
t42 = g(1) * t5 + g(2) * t6;
t36 = t5 * pkin(3);
t35 = t6 * pkin(3);
t28 = t32 * pkin(1) + t30 * qJ(2);
t27 = t32 * pkin(2) + t28;
t26 = -g(1) * t35 + g(2) * t36;
t25 = -t30 * pkin(1) + t32 * qJ(2);
t24 = -t6 * rSges(4,1) + t5 * rSges(4,2);
t23 = t5 * rSges(4,1) + t6 * rSges(4,2);
t20 = -t30 * pkin(2) + t25;
t18 = g(1) * (t35 + t20) + g(2) * (-t36 + t27);
t1 = [-m(2) * (g(1) * (-t30 * rSges(2,1) - t32 * rSges(2,2)) + g(2) * (t32 * rSges(2,1) - t30 * rSges(2,2))) - m(3) * (g(1) * (-t30 * rSges(3,1) + t32 * rSges(3,3) + t25) + g(2) * (t32 * rSges(3,1) + t30 * rSges(3,3) + t28)) - m(4) * (g(1) * (t20 - t24) + g(2) * (-t23 + t27)) - m(5) * (t18 + t46) - m(6) * (t18 + t49), (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t30 - g(2) * t32), -m(4) * (g(1) * t24 + g(2) * t23) - m(5) * (t26 - t46) - m(6) * (t26 - t49), (m(5) * t22 + m(6) * t19) * g(3) + t42 * (-m(5) * (rSges(5,1) * t15 + rSges(5,2) * t16) - m(6) * (t43 * t15 - t41 * t16)), -m(6) * (g(3) * t16 - t42 * t15)];
taug = t1(:);

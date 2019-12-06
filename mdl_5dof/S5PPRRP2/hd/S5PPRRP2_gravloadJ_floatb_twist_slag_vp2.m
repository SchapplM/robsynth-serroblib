% Calculate Gravitation load on the joints for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:32
% EndTime: 2019-12-05 15:08:33
% DurationCPUTime: 0.26s
% Computational Cost: add. (148->35), mult. (206->47), div. (0->0), fcn. (180->6), ass. (0->24)
t41 = mrSges(5,1) + mrSges(6,1);
t40 = -mrSges(5,2) + mrSges(6,3);
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t39 = t40 * t16 + t41 * t17 + mrSges(4,1);
t38 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t36 = -m(5) - m(6);
t33 = m(6) * pkin(4) + t41;
t32 = -m(6) * qJ(5) - t40;
t13 = pkin(8) + qJ(3);
t11 = sin(t13);
t29 = g(3) * t11;
t14 = sin(pkin(7));
t28 = t14 * t16;
t27 = t14 * t17;
t15 = cos(pkin(7));
t26 = t15 * t16;
t25 = t15 * t17;
t24 = m(3) + m(4) - t36;
t20 = pkin(4) * t17 + qJ(5) * t16;
t12 = cos(t13);
t3 = t12 * t26 - t27;
t1 = t12 * t28 + t25;
t2 = [(-m(2) - t24) * g(3), (-g(1) * t14 + g(2) * t15) * t24, (t36 * (t12 * pkin(3) + t11 * pkin(6)) + (-m(6) * t20 - t39) * t12 + t38 * t11) * g(3) + (g(1) * t15 + g(2) * t14) * ((-m(6) * (-pkin(3) - t20) + m(5) * pkin(3) + t39) * t11 + (t36 * pkin(6) + t38) * t12), (t33 * t16 + t32 * t17) * t29 + (t32 * (t12 * t27 - t26) + t33 * t1) * g(2) + (t32 * (t12 * t25 + t28) + t33 * t3) * g(1), (-g(1) * t3 - g(2) * t1 - t16 * t29) * m(6)];
taug = t2(:);

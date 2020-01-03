% Calculate Gravitation load on the joints for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:14
% DurationCPUTime: 0.25s
% Computational Cost: add. (77->33), mult. (189->46), div. (0->0), fcn. (172->6), ass. (0->20)
t37 = mrSges(4,1) + mrSges(5,1);
t36 = -mrSges(4,2) + mrSges(5,3);
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t35 = t36 * t13 + t37 * t15 + mrSges(3,1);
t34 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t32 = -m(4) - m(5);
t29 = m(5) * pkin(3) + t37;
t28 = -m(5) * qJ(4) - t36;
t14 = sin(qJ(2));
t26 = g(3) * t14;
t16 = cos(qJ(2));
t24 = t16 * t13;
t23 = t16 * t15;
t19 = pkin(3) * t15 + qJ(4) * t13;
t12 = cos(pkin(6));
t11 = sin(pkin(6));
t3 = -t11 * t15 + t12 * t24;
t1 = t11 * t24 + t12 * t15;
t2 = [(-m(2) - m(3) + t32) * g(3), (t32 * (t16 * pkin(2) + t14 * pkin(5)) + (-m(5) * t19 - t35) * t16 + t34 * t14) * g(3) + (g(1) * t12 + g(2) * t11) * ((-m(5) * (-pkin(2) - t19) + m(4) * pkin(2) + t35) * t14 + (t32 * pkin(5) + t34) * t16), (t29 * t13 + t28 * t15) * t26 + (t28 * (t11 * t23 - t12 * t13) + t29 * t1) * g(2) + (t28 * (t11 * t13 + t12 * t23) + t29 * t3) * g(1), (-g(1) * t3 - g(2) * t1 - t13 * t26) * m(5)];
taug = t2(:);

% Calculate Gravitation load on the joints for
% S4PRRP3
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:44
% DurationCPUTime: 0.13s
% Computational Cost: add. (81->22), mult. (86->22), div. (0->0), fcn. (61->4), ass. (0->12)
t11 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t14 = mrSges(4,2) + mrSges(5,2);
t6 = sin(qJ(3));
t7 = cos(qJ(3));
t20 = -t11 * t7 + t14 * t6;
t19 = m(4) + m(5);
t16 = t19 * pkin(2) + mrSges(3,1) - t20;
t15 = mrSges(3,2) + m(5) * (-qJ(4) - pkin(5)) - mrSges(5,3) - m(4) * pkin(5) - mrSges(4,3);
t4 = pkin(6) + qJ(2);
t3 = cos(t4);
t2 = sin(t4);
t1 = [(-m(2) - m(3) - t19) * g(3), (t15 * t2 - t16 * t3) * g(2) + (t15 * t3 + t16 * t2) * g(1), t20 * g(3) + (g(1) * t3 + g(2) * t2) * (t11 * t6 + t14 * t7), (-g(1) * t2 + g(2) * t3) * m(5)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:09
% EndTime: 2019-12-31 17:18:11
% DurationCPUTime: 0.38s
% Computational Cost: add. (122->54), mult. (261->68), div. (0->0), fcn. (236->6), ass. (0->31)
t43 = -mrSges(4,1) - mrSges(5,1);
t41 = mrSges(4,2) + mrSges(5,2);
t47 = mrSges(4,3) + mrSges(5,3);
t10 = sin(qJ(3));
t13 = cos(qJ(3));
t5 = pkin(3) * t13 + pkin(2);
t46 = -m(4) * pkin(2) - m(5) * t5 + t41 * t10 + t43 * t13;
t9 = -qJ(4) - pkin(6);
t45 = -m(4) * pkin(6) + m(5) * t9 - t47;
t12 = sin(qJ(1));
t15 = cos(qJ(1));
t44 = g(1) * t15 + g(2) * t12;
t36 = m(5) * pkin(3);
t42 = mrSges(2,2) - mrSges(3,3);
t40 = -m(3) - m(4) - m(5);
t39 = t36 - t43;
t11 = sin(qJ(2));
t14 = cos(qJ(2));
t21 = mrSges(3,1) * t14 - mrSges(3,2) * t11;
t37 = t47 * t11 + mrSges(2,1) + t21;
t31 = t10 * t15;
t28 = t12 * t10;
t27 = t12 * t14;
t26 = t14 * t15;
t23 = pkin(2) * t14 + pkin(6) * t11;
t22 = -t11 * t9 + t14 * t5;
t3 = -t10 * t26 + t12 * t13;
t1 = t10 * t27 + t13 * t15;
t4 = t13 * t26 + t28;
t2 = -t13 * t27 + t31;
t6 = [(-t28 * t36 + t43 * t4 + t40 * (t15 * pkin(1) + t12 * pkin(5)) - t41 * t3 + t42 * t12 + (-m(4) * t23 - m(5) * t22 - t37) * t15) * g(2) + (-t31 * t36 + t43 * t2 - t41 * t1 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t23) - m(5) * (-pkin(1) - t22) + t37) * t12 + (t40 * pkin(5) + t42) * t15) * g(1), -g(3) * t21 + (t46 * g(3) + t44 * (mrSges(3,2) + t45)) * t14 + (t45 * g(3) + t44 * (mrSges(3,1) - t46)) * t11, (t10 * t39 + t13 * t41) * g(3) * t11 + (t1 * t39 - t2 * t41) * g(2) + (-t3 * t39 + t4 * t41) * g(1), (g(3) * t14 - t44 * t11) * m(5)];
taug = t6(:);

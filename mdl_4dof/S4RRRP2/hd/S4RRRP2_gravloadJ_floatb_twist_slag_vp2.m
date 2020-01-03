% Calculate Gravitation load on the joints for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:12:54
% DurationCPUTime: 0.18s
% Computational Cost: add. (139->38), mult. (140->41), div. (0->0), fcn. (105->6), ass. (0->21)
t22 = cos(qJ(3));
t20 = sin(qJ(3));
t31 = mrSges(4,2) + mrSges(5,2);
t28 = t31 * t20;
t32 = mrSges(4,1) + mrSges(5,1);
t37 = -t32 * t22 - mrSges(3,1) + t28;
t35 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t18 = qJ(1) + qJ(2);
t15 = sin(t18);
t16 = cos(t18);
t30 = t16 * pkin(2) + t15 * pkin(6);
t14 = pkin(3) * t22 + pkin(2);
t19 = -qJ(4) - pkin(6);
t29 = t16 * t14 - t15 * t19;
t27 = m(5) * pkin(3) + t32;
t25 = t35 * t15 + t37 * t16;
t24 = (-m(4) * pkin(6) + m(5) * t19 + t35) * t16 + (m(4) * pkin(2) + m(5) * t14 - t37) * t15;
t23 = cos(qJ(1));
t21 = sin(qJ(1));
t17 = t23 * pkin(1);
t1 = [(t21 * mrSges(2,2) - m(4) * (t17 + t30) - m(5) * (t17 + t29) + (-m(3) * pkin(1) - mrSges(2,1)) * t23 + t25) * g(2) + (mrSges(2,2) * t23 + (mrSges(2,1) + (m(3) + m(4) + m(5)) * pkin(1)) * t21 + t24) * g(1), (-m(4) * t30 - m(5) * t29 + t25) * g(2) + t24 * g(1), (-t27 * t22 + t28) * g(3) + (g(1) * t16 + g(2) * t15) * (t27 * t20 + t31 * t22), (-g(1) * t15 + g(2) * t16) * m(5)];
taug = t1(:);

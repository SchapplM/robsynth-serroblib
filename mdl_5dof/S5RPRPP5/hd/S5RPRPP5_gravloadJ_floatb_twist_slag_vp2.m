% Calculate Gravitation load on the joints for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:09
% DurationCPUTime: 0.35s
% Computational Cost: add. (104->57), mult. (206->58), div. (0->0), fcn. (157->4), ass. (0->22)
t44 = mrSges(4,1) + mrSges(5,1);
t10 = sin(qJ(3));
t12 = cos(qJ(3));
t43 = -t44 * t12 + (mrSges(4,2) - mrSges(5,3)) * t10;
t11 = sin(qJ(1));
t28 = g(1) * t11;
t42 = pkin(3) * t12 + qJ(4) * t10;
t41 = (mrSges(4,2) - mrSges(6,2)) * t12 + t44 * t10;
t40 = -m(3) - m(4);
t39 = m(5) + m(6);
t13 = cos(qJ(1));
t27 = g(2) * t13;
t38 = -t28 + t27;
t6 = t12 * qJ(4);
t37 = mrSges(2,2) - mrSges(3,3) - (-m(5) * qJ(4) - mrSges(5,3)) * t12 - m(6) * (pkin(4) * t10 - t6) - t10 * mrSges(6,1) - t41;
t36 = m(6) * qJ(5) - mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t31 = t13 * pkin(1) + t11 * qJ(2);
t30 = pkin(3) * t10;
t23 = t13 * pkin(6) + t31;
t20 = m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1);
t7 = t13 * qJ(2);
t1 = [(-m(3) * t31 - m(4) * t23 - t39 * (t11 * t30 + t23) + t36 * t13 + t37 * t11) * g(2) + (t40 * t7 - t39 * (t13 * t30 + t7) + t37 * t13 + (m(3) * pkin(1) + (-m(4) - t39) * (-pkin(1) - pkin(6)) - t36) * t11) * g(1), t38 * (t39 - t40), (m(5) * t42 - (-m(6) * qJ(4) - mrSges(6,2)) * t10 - t20 * t12 - t43) * t27 + (-m(5) * (t6 - t30) - t12 * mrSges(5,3) - m(6) * t6 - t20 * t10 + t41) * g(3) + (-t39 * t42 - t10 * mrSges(6,2) - (m(6) * pkin(4) + mrSges(6,1)) * t12 + t43) * t28, (-g(3) * t10 - t12 * t38) * t39, (g(1) * t13 + g(2) * t11) * m(6)];
taug = t1(:);

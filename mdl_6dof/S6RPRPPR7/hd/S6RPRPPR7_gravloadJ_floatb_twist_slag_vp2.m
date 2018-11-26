% Calculate Gravitation load on the joints for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:55:26
% EndTime: 2018-11-23 15:55:27
% DurationCPUTime: 0.53s
% Computational Cost: add. (237->88), mult. (325->96), div. (0->0), fcn. (265->8), ass. (0->44)
t71 = -mrSges(5,1) + mrSges(6,2);
t19 = qJ(3) + pkin(9);
t14 = sin(t19);
t15 = cos(t19);
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t59 = pkin(3) * t25;
t70 = -m(5) * t59 - mrSges(4,1) * t25 + mrSges(4,2) * t22 + t71 * t15 + (mrSges(5,2) - mrSges(6,3)) * t14;
t23 = sin(qJ(1));
t57 = g(1) * t23;
t69 = -t22 * mrSges(4,1) - t25 * mrSges(4,2) - t15 * mrSges(5,2) + t71 * t14;
t68 = pkin(4) * t15 + qJ(5) * t14 + t59;
t67 = -m(3) - m(4);
t66 = m(6) + m(7);
t26 = cos(qJ(1));
t56 = g(2) * t26;
t63 = -t57 + t56;
t12 = t15 * qJ(5);
t62 = mrSges(2,2) - mrSges(3,3) - (-m(6) * qJ(5) - mrSges(6,3)) * t15 - m(7) * (t14 * pkin(8) - t12) - t14 * mrSges(7,3) + t69;
t61 = m(7) * pkin(5) + mrSges(2,1) + mrSges(6,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t55 = g(3) * t14;
t54 = t14 * pkin(4);
t53 = t22 * pkin(3);
t21 = sin(qJ(6));
t51 = t23 * t21;
t24 = cos(qJ(6));
t50 = t23 * t24;
t49 = t26 * t21;
t48 = t26 * t24;
t47 = t26 * pkin(1) + t23 * qJ(2);
t45 = -m(5) - t66;
t17 = t26 * qJ(2);
t43 = -t23 * pkin(1) + t17;
t42 = t12 - t53;
t40 = m(7) * (-pkin(4) - pkin(8)) - mrSges(7,3);
t35 = t21 * mrSges(7,1) + t24 * mrSges(7,2);
t20 = -qJ(4) - pkin(7);
t33 = t23 * t20 + t26 * t53 + t43;
t32 = -t26 * t20 + t23 * t53 + t47;
t4 = -t15 * t51 + t48;
t3 = -t15 * t50 - t49;
t2 = -t15 * t49 - t50;
t1 = -t15 * t48 + t51;
t5 = [(-m(5) * t32 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + t67 * t47 - t66 * (t23 * t54 + t32) + (-m(4) * pkin(7) - t61) * t26 + t62 * t23) * g(2) + (-m(3) * t43 - m(4) * t17 - m(5) * t33 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t66 * (t26 * t54 + t33) + (-m(4) * (-pkin(1) - pkin(7)) + t61) * t23 + t62 * t26) * g(1), t63 * (-t45 - t67) (m(6) * t68 + m(7) * t59 - t40 * t15 - (-m(7) * qJ(5) - t35) * t14 - t70) * t56 + (m(5) * t53 - m(6) * (t42 - t54) - m(7) * t42 - t40 * t14 + (-mrSges(6,3) - t35) * t15 - t69) * g(3) + (-t66 * t68 - (m(7) * pkin(8) + mrSges(7,3)) * t15 - t35 * t14 + t70) * t57 (g(1) * t26 + g(2) * t23) * t45 (-t15 * t63 - t55) * t66, -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - (mrSges(7,1) * t24 - mrSges(7,2) * t21) * t55];
taug  = t5(:);

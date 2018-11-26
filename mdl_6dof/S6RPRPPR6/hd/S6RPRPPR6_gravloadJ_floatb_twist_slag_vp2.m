% Calculate Gravitation load on the joints for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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

function taug = S6RPRPPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:54:42
% EndTime: 2018-11-23 15:54:42
% DurationCPUTime: 0.56s
% Computational Cost: add. (273->86), mult. (344->91), div. (0->0), fcn. (288->10), ass. (0->43)
t65 = m(6) + m(7);
t73 = m(5) + t65;
t72 = mrSges(5,2) - mrSges(7,3);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t60 = -g(1) * t24 + g(2) * t26;
t18 = qJ(3) + pkin(9);
t11 = sin(t18);
t13 = cos(t18);
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t69 = t23 * mrSges(4,1) + t11 * mrSges(5,1) + t25 * mrSges(4,2) + t72 * t13;
t67 = -m(3) - m(4);
t66 = -m(5) - m(7);
t21 = -qJ(4) - pkin(7);
t51 = t23 * pkin(3);
t64 = t24 * t21 + t26 * t51;
t17 = pkin(10) + qJ(6);
t10 = sin(t17);
t12 = cos(t17);
t19 = sin(pkin(10));
t20 = cos(pkin(10));
t31 = m(6) * pkin(4) + t20 * mrSges(6,1) - t19 * mrSges(6,2);
t8 = t20 * pkin(5) + pkin(4);
t63 = m(7) * t8 + t12 * mrSges(7,1) - t10 * mrSges(7,2) + t31;
t40 = m(6) * qJ(5) + mrSges(6,3);
t22 = -pkin(8) - qJ(5);
t49 = t13 * t22;
t58 = -t31 * t11 + t40 * t13 + mrSges(2,2) - mrSges(3,3) - m(7) * (t11 * t8 + t49) - t69;
t57 = t20 * mrSges(6,2) + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t19;
t48 = t24 * t10;
t47 = t24 * t12;
t46 = t26 * t10;
t45 = t26 * t12;
t44 = t26 * pkin(1) + t24 * qJ(2);
t42 = t24 * t51 + t44;
t15 = t26 * qJ(2);
t41 = -t24 * pkin(1) + t15;
t4 = t11 * t45 - t48;
t3 = t11 * t46 + t47;
t2 = t11 * t47 + t46;
t1 = -t11 * t48 + t45;
t5 = [(-m(6) * t42 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t67 * t44 + t66 * (-t26 * t21 + t42) + (-m(4) * pkin(7) + m(6) * t21 - t57) * t26 + t58 * t24) * g(2) + (-m(3) * t41 - m(4) * t15 - m(6) * (t15 + t64) - t4 * mrSges(7,1) + t3 * mrSges(7,2) + t66 * (t41 + t64) + (-m(4) * (-pkin(1) - pkin(7)) + m(6) * pkin(1) + t57) * t24 + t58 * t26) * g(1), t60 * (t73 - t67) (m(5) * t51 - m(6) * (t13 * qJ(5) - t51) - t13 * mrSges(6,3) - m(7) * (-t49 - t51) + t63 * t11 + t69) * g(3) + t60 * (-mrSges(4,2) * t23 + (mrSges(5,1) + t63) * t13 + (-m(7) * t22 + t40 - t72) * t11 + (t73 * pkin(3) + mrSges(4,1)) * t25) -(g(1) * t26 + g(2) * t24) * t73 (-g(3) * t11 - t13 * t60) * t65, -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t10 - mrSges(7,2) * t12) * t13];
taug  = t5(:);

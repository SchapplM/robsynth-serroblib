% Calculate Gravitation load on the joints for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:53
% EndTime: 2019-03-09 08:13:55
% DurationCPUTime: 0.72s
% Computational Cost: add. (243->100), mult. (448->100), div. (0->0), fcn. (384->8), ass. (0->47)
t84 = mrSges(3,1) + mrSges(4,1);
t83 = mrSges(5,2) - mrSges(7,3);
t82 = mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t74 = m(6) + m(7);
t51 = m(5) + t74;
t81 = m(4) + t51;
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t72 = g(1) * t26 + g(2) * t24;
t22 = -pkin(8) - qJ(5);
t66 = -pkin(2) - pkin(3);
t80 = -m(5) * t66 - m(6) * (-qJ(5) + t66) + mrSges(6,3) - m(7) * (t22 + t66) - t83;
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t77 = t82 * t23 + t84 * t25;
t75 = -m(5) - m(7);
t20 = sin(pkin(9));
t21 = cos(pkin(9));
t38 = -t21 * mrSges(6,1) + t20 * mrSges(6,2);
t34 = m(6) * pkin(4) - t38;
t69 = -t34 * t23 + t83 * t25 - t77;
t68 = m(6) * qJ(4) + t21 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t20;
t62 = g(3) * t25;
t16 = t25 * pkin(2);
t61 = t22 * t25;
t59 = t23 * t26;
t19 = pkin(9) + qJ(6);
t11 = sin(t19);
t58 = t24 * t11;
t12 = cos(t19);
t57 = t24 * t12;
t55 = t25 * t26;
t13 = t23 * qJ(3);
t54 = t16 + t13;
t53 = t26 * pkin(1) + t24 * pkin(7);
t50 = t25 * pkin(3) + t54;
t48 = pkin(2) * t55 + t26 * t13 + t53;
t47 = -pkin(1) - t13;
t45 = pkin(3) * t55 + t48;
t10 = pkin(5) * t21 + pkin(4);
t33 = m(7) * t10 + t12 * mrSges(7,1) - t11 * mrSges(7,2);
t17 = t26 * pkin(7);
t4 = t12 * t59 - t58;
t3 = -t11 * t59 - t57;
t2 = -t11 * t26 - t23 * t57;
t1 = -t12 * t26 + t23 * t58;
t5 = [(-m(3) * t53 - m(4) * t48 - m(6) * t45 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + t75 * (-t24 * qJ(4) + t45) + t68 * t24 + (-mrSges(2,1) - (m(6) * qJ(5) + mrSges(6,3)) * t25 - m(7) * (t10 * t23 - t61) + t69) * t26) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) + t75 * (-qJ(4) * t26 + t17) + (-m(3) - m(4) - m(6)) * t17 + t68 * t26 + (mrSges(2,1) - m(4) * (t47 - t16) - m(5) * t47 + (-m(6) * (-pkin(4) - qJ(3)) - t38 - m(7) * (-qJ(3) - t10)) * t23 + (m(3) + t74) * pkin(1) + t80 * t25 + t77) * t24) * g(1) (-m(4) * t54 - m(5) * t50 - m(6) * (qJ(5) * t25 + t50) - t25 * mrSges(6,3) - m(7) * (t50 - t61) - t33 * t23 + t69) * g(3) + ((m(4) * pkin(2) + t80 + t84) * t23 + (-qJ(3) * t81 - t33 - t34 - t82) * t25) * t72 (-t72 * t23 + t62) * t81 (g(1) * t24 - g(2) * t26) * t51 (-t23 * g(3) - t72 * t25) * t74, -g(1) * (mrSges(7,1) * t3 - mrSges(7,2) * t4) - g(2) * (-mrSges(7,1) * t1 + mrSges(7,2) * t2) - (mrSges(7,1) * t11 + mrSges(7,2) * t12) * t62];
taug  = t5(:);

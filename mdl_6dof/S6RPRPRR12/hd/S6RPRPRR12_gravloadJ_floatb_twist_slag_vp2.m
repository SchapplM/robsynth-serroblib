% Calculate Gravitation load on the joints for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:18
% EndTime: 2019-03-09 04:18:20
% DurationCPUTime: 0.74s
% Computational Cost: add. (236->104), mult. (421->121), div. (0->0), fcn. (367->8), ass. (0->53)
t33 = -pkin(9) - pkin(8);
t77 = mrSges(6,3) + mrSges(7,3);
t82 = -m(6) * pkin(8) + m(7) * t33 - t77;
t80 = -m(7) * pkin(5) - mrSges(6,1);
t28 = sin(qJ(3));
t31 = cos(qJ(3));
t79 = t31 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * t28;
t78 = -m(3) - m(4);
t26 = qJ(5) + qJ(6);
t19 = sin(t26);
t20 = cos(t26);
t27 = sin(qJ(5));
t30 = cos(qJ(5));
t76 = t27 * mrSges(6,1) + t19 * mrSges(7,1) + t30 * mrSges(6,2) + t20 * mrSges(7,2);
t49 = m(5) + m(6) + m(7);
t22 = t31 * qJ(4);
t66 = pkin(5) * t27;
t75 = mrSges(2,2) - mrSges(3,3) - (-m(5) * qJ(4) - mrSges(5,3)) * t31 + m(6) * t22 - m(7) * (-t31 * t66 - t22) - t79 + t82 * t28;
t61 = t30 * pkin(5);
t74 = m(6) * pkin(4) + m(7) * (pkin(4) + t61) + mrSges(2,1) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t73 = -m(7) * t66 - t76;
t72 = -m(6) * (-pkin(3) - pkin(8)) - m(7) * (-pkin(3) + t33) + t77;
t32 = cos(qJ(1));
t64 = g(2) * t32;
t29 = sin(qJ(1));
t65 = g(1) * t29;
t71 = -t65 + t64;
t56 = t32 * t20;
t5 = t19 * t29 - t31 * t56;
t57 = t32 * t19;
t6 = -t20 * t29 - t31 * t57;
t68 = -mrSges(7,1) * t5 + mrSges(7,2) * t6;
t58 = t29 * t31;
t7 = -t20 * t58 - t57;
t8 = -t19 * t58 + t56;
t67 = mrSges(7,1) * t7 - mrSges(7,2) * t8;
t63 = g(3) * t28;
t62 = t28 * pkin(3);
t60 = mrSges(7,2) * t19;
t55 = t32 * t27;
t54 = t32 * t30;
t51 = pkin(1) * t32 + qJ(2) * t29;
t50 = qJ(4) * t28;
t48 = pkin(7) * t32 + t51;
t43 = mrSges(4,1) * t31 - mrSges(4,2) * t28;
t39 = -t31 * mrSges(5,2) + t28 * mrSges(5,3);
t38 = -t27 * t29 + t31 * t54;
t11 = -t30 * t58 - t55;
t23 = t32 * qJ(2);
t13 = t28 * t20 * mrSges(7,1);
t12 = -t27 * t58 + t54;
t10 = -t29 * t30 - t31 * t55;
t1 = [(-m(3) * t51 - m(4) * t48 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t49 * (t29 * t62 + t48) - t74 * t32 + t75 * t29) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) + t38 * mrSges(6,2) - t5 * mrSges(7,2) + t78 * t23 - t49 * (t32 * t62 + t23) + (m(3) * pkin(1) + (-m(4) - t49) * (-pkin(1) - pkin(7)) + t74) * t29 + t75 * t32) * g(1), t71 * (t49 - t78) -t43 * t65 + (-t49 * (pkin(3) * t58 + t29 * t50) + (t73 * t28 + t31 * t82 - t39) * t29) * g(1) + (m(5) * t50 + t39 + t43 + (m(5) * pkin(3) + t72) * t31 + (m(6) * qJ(4) - m(7) * (-qJ(4) - t66) + t76) * t28) * t64 + (m(5) * t62 + t72 * t28 - t49 * t22 + (-mrSges(5,3) + t73) * t31 + t79) * g(3) (-t31 * t71 - t63) * t49 -(mrSges(6,1) * t30 - mrSges(6,2) * t27) * t63 - g(3) * (t13 + (m(7) * t61 - t60) * t28) + (-t10 * mrSges(6,2) + t38 * t80 - t68) * g(2) + (t12 * mrSges(6,2) + t11 * t80 - t67) * g(1), -g(1) * t67 - g(2) * t68 - g(3) * (-t28 * t60 + t13)];
taug  = t1(:);

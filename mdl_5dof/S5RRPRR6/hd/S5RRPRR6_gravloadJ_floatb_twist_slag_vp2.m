% Calculate Gravitation load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:50
% EndTime: 2022-01-20 11:16:52
% DurationCPUTime: 0.40s
% Computational Cost: add. (334->84), mult. (316->99), div. (0->0), fcn. (294->10), ass. (0->48)
t90 = -mrSges(4,3) + mrSges(3,2);
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t89 = -t46 * mrSges(4,1) - mrSges(3,1) + (mrSges(4,2) - mrSges(5,3) - mrSges(6,3)) * t45;
t87 = m(6) * pkin(4) + mrSges(5,1);
t44 = qJ(1) + qJ(2);
t39 = sin(t44);
t70 = t45 * (-pkin(8) - pkin(7));
t41 = cos(t44);
t47 = sin(qJ(4));
t74 = t41 * t47;
t86 = pkin(4) * t74 + t39 * t70;
t43 = qJ(4) + qJ(5);
t38 = sin(t43);
t40 = cos(t43);
t77 = t39 * t46;
t10 = t41 * t38 - t40 * t77;
t49 = cos(qJ(4));
t68 = t46 * t47;
t17 = t39 * t68 + t41 * t49;
t67 = t46 * t49;
t18 = -t39 * t67 + t74;
t37 = t49 * pkin(4) + pkin(3);
t80 = pkin(7) * t45;
t9 = t38 * t77 + t41 * t40;
t85 = -t18 * mrSges(5,1) - t10 * mrSges(6,1) - t17 * mrSges(5,2) - t9 * mrSges(6,2) + t90 * t41 + (-m(5) * (-pkin(3) * t46 - pkin(2) - t80) - m(6) * (-t37 * t46 - pkin(2)) - t89) * t39;
t75 = t41 * t46;
t11 = -t38 * t75 + t39 * t40;
t12 = t39 * t38 + t40 * t75;
t19 = t39 * t49 - t41 * t68;
t76 = t39 * t47;
t20 = t41 * t67 + t76;
t84 = -t20 * mrSges(5,1) - t12 * mrSges(6,1) - t19 * mrSges(5,2) - t11 * mrSges(6,2) + t39 * t90 + t89 * t41;
t82 = -t9 * mrSges(6,1) + t10 * mrSges(6,2);
t81 = t11 * mrSges(6,1) - t12 * mrSges(6,2);
t79 = g(3) * t45;
t48 = sin(qJ(1));
t78 = t48 * pkin(1);
t50 = cos(qJ(1));
t42 = t50 * pkin(1);
t66 = t41 * pkin(2) + t39 * qJ(3);
t32 = t41 * qJ(3);
t62 = t32 - t78;
t61 = -t39 * pkin(2) + t32;
t59 = pkin(3) * t75 + t41 * t80 + t66;
t57 = -mrSges(6,1) * t38 - mrSges(6,2) * t40;
t54 = pkin(4) * t76 + t37 * t75 - t41 * t70 + t66;
t1 = [(-t50 * mrSges(2,1) + t48 * mrSges(2,2) - m(3) * t42 - m(4) * (t42 + t66) - m(5) * (t42 + t59) - m(6) * (t42 + t54) + t84) * g(2) + (t48 * mrSges(2,1) + t50 * mrSges(2,2) + m(3) * t78 - m(4) * (t61 - t78) - m(5) * t62 - m(6) * (t62 + t86) + t85) * g(1), (-m(4) * t66 - m(5) * t59 - m(6) * t54 + t84) * g(2) + (-m(4) * t61 - m(5) * t32 - m(6) * (t32 + t86) + t85) * g(1), (-g(1) * t39 + g(2) * t41) * (m(4) + m(5) + m(6)), (mrSges(5,2) * t49 + t87 * t47 - t57) * t79 + (-t18 * mrSges(5,2) + t87 * t17 - t82) * g(2) + (t20 * mrSges(5,2) - t87 * t19 - t81) * g(1), -g(1) * t81 - g(2) * t82 - t57 * t79];
taug = t1(:);

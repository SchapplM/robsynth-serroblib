% Calculate Gravitation load on the joints for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:34
% EndTime: 2019-03-08 19:21:35
% DurationCPUTime: 0.59s
% Computational Cost: add. (401->83), mult. (1003->125), div. (0->0), fcn. (1203->12), ass. (0->47)
t83 = m(6) + m(7);
t45 = sin(qJ(6));
t48 = cos(qJ(6));
t80 = -m(7) * pkin(5) - mrSges(7,1) * t48 + mrSges(7,2) * t45 - mrSges(6,1);
t78 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t81 = mrSges(3,1) + mrSges(4,1);
t79 = mrSges(3,2) - mrSges(4,3);
t69 = m(5) + t83;
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t77 = t78 * t46 + t80 * t49 - mrSges(5,1);
t76 = t45 * mrSges(7,1) + t48 * mrSges(7,2) + t83 * pkin(8) - mrSges(5,2) + mrSges(6,3);
t42 = sin(pkin(6));
t75 = t42 * t46;
t47 = sin(qJ(2));
t74 = t42 * t47;
t73 = t42 * t49;
t50 = cos(qJ(2));
t72 = t42 * t50;
t71 = pkin(2) * t72 + qJ(3) * t74;
t70 = cos(pkin(6));
t68 = pkin(3) * t72 + t71;
t67 = m(4) + t69;
t41 = sin(pkin(10));
t44 = cos(pkin(10));
t62 = t50 * t70;
t28 = t41 * t47 - t44 * t62;
t63 = t47 * t70;
t29 = t41 * t50 + t44 * t63;
t40 = sin(pkin(11));
t43 = cos(pkin(11));
t65 = -t28 * t43 + t29 * t40;
t10 = t28 * t40 + t29 * t43;
t30 = t41 * t62 + t44 * t47;
t31 = -t41 * t63 + t44 * t50;
t64 = -t30 * t43 + t31 * t40;
t14 = t30 * t40 + t31 * t43;
t61 = -t28 * pkin(2) + t29 * qJ(3);
t60 = -t30 * pkin(2) + t31 * qJ(3);
t58 = -t28 * pkin(3) + t61;
t57 = -t30 * pkin(3) + t60;
t23 = -t40 * t72 + t43 * t74;
t22 = (t40 * t47 + t43 * t50) * t42;
t16 = t23 * t49 - t70 * t46;
t4 = t14 * t49 - t41 * t75;
t2 = t10 * t49 + t44 * t75;
t1 = [(-m(2) - m(3) - t67) * g(3) (-m(4) * t71 - m(5) * t68 - t83 * (t22 * pkin(4) + t68) + (t79 * t47 - t81 * t50) * t42 + t76 * t23 + t77 * t22) * g(3) + (-m(4) * t61 - m(5) * t58 - t83 * (pkin(4) * t65 + t58) + t79 * t29 + t81 * t28 + t77 * t65 + t76 * t10) * g(2) + (-m(4) * t60 - m(5) * t57 - t83 * (pkin(4) * t64 + t57) + t79 * t31 + t81 * t30 + t77 * t64 + t76 * t14) * g(1) (-g(1) * t30 - g(2) * t28 + g(3) * t72) * t67 ((g(1) * t41 - g(2) * t44) * t42 + g(3) * t70) * t69 (t78 * t16 + t80 * (-t23 * t46 - t70 * t49)) * g(3) + (t78 * t2 + t80 * (-t10 * t46 + t44 * t73)) * g(2) + (t78 * t4 + t80 * (-t14 * t46 - t41 * t73)) * g(1), -g(1) * ((-t4 * t45 + t48 * t64) * mrSges(7,1) + (-t4 * t48 - t45 * t64) * mrSges(7,2)) - g(2) * ((-t2 * t45 + t48 * t65) * mrSges(7,1) + (-t2 * t48 - t45 * t65) * mrSges(7,2)) - g(3) * ((-t16 * t45 + t22 * t48) * mrSges(7,1) + (-t16 * t48 - t22 * t45) * mrSges(7,2))];
taug  = t1(:);

% Calculate Gravitation load on the joints for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:40
% EndTime: 2019-03-08 21:22:43
% DurationCPUTime: 1.02s
% Computational Cost: add. (553->98), mult. (990->145), div. (0->0), fcn. (1123->12), ass. (0->51)
t46 = cos(qJ(5));
t107 = -m(6) * pkin(4) - m(7) * (pkin(5) * t46 + pkin(4)) - mrSges(5,1);
t90 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t106 = mrSges(6,1) + mrSges(7,1);
t97 = -mrSges(6,2) - mrSges(7,2);
t92 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t89 = m(7) * pkin(5);
t93 = -t89 - t106;
t39 = qJ(3) + pkin(11);
t37 = sin(t39);
t38 = cos(t39);
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t105 = -m(4) * pkin(2) - t47 * mrSges(4,1) + t44 * mrSges(4,2) + t107 * t38 + t90 * t37 - mrSges(3,1);
t95 = m(5) + m(6) + m(7);
t45 = sin(qJ(2));
t48 = cos(qJ(2));
t71 = sin(pkin(10));
t73 = cos(pkin(6));
t57 = t73 * t71;
t72 = cos(pkin(10));
t26 = -t45 * t57 + t72 * t48;
t40 = sin(pkin(6));
t67 = t40 * t71;
t104 = -t26 * t44 + t47 * t67;
t79 = t40 * t45;
t103 = -t44 * t79 + t73 * t47;
t58 = t73 * t72;
t24 = t45 * t58 + t71 * t48;
t68 = t40 * t72;
t96 = -t24 * t44 - t47 * t68;
t43 = sin(qJ(5));
t94 = t106 * t46 + t97 * t43 - t107;
t87 = t24 * t43;
t85 = t26 * t43;
t81 = t38 * t43;
t80 = t38 * t46;
t78 = t40 * t48;
t77 = t43 * t48;
t76 = t46 * t48;
t42 = -qJ(4) - pkin(8);
t36 = pkin(3) * t47 + pkin(2);
t25 = t72 * t45 + t48 * t57;
t23 = t71 * t45 - t48 * t58;
t20 = t73 * t37 + t38 * t79;
t19 = t37 * t79 - t73 * t38;
t12 = t26 * t38 + t37 * t67;
t11 = t26 * t37 - t38 * t67;
t10 = t24 * t38 - t37 * t68;
t9 = t24 * t37 + t38 * t68;
t1 = [(-m(2) - m(3) - m(4) - t95) * g(3) (-t87 * t89 - t95 * (-t23 * t36 - t24 * t42) - t106 * (-t23 * t80 + t87) + t97 * (t23 * t81 + t24 * t46) + t92 * t24 - t105 * t23) * g(2) + (-t85 * t89 - t106 * (-t25 * t80 + t85) - t95 * (-t25 * t36 - t26 * t42) + t97 * (t25 * t81 + t26 * t46) + t92 * t26 - t105 * t25) * g(1) + (-t95 * t36 * t78 + (t105 * t48 + (-t106 * t76 - t97 * t77) * t38 + (t95 * t42 + t93 * t43 + t97 * t46 + t92) * t45) * t40) * g(3) (-t103 * mrSges(4,1) - (-t73 * t44 - t47 * t79) * mrSges(4,2) + t90 * t20 + t94 * t19) * g(3) + (-(-t24 * t47 + t44 * t68) * mrSges(4,2) - mrSges(4,1) * t96 + t94 * t9 + t90 * t10) * g(2) + (-t104 * mrSges(4,1) - (-t26 * t47 - t44 * t67) * mrSges(4,2) + t90 * t12 + t94 * t11) * g(1) + (-g(1) * t104 - t96 * g(2) - g(3) * t103) * t95 * pkin(3), t95 * (-g(1) * t25 - g(2) * t23 + g(3) * t78) (t97 * (-t20 * t46 + t40 * t77) + t93 * (-t20 * t43 - t40 * t76)) * g(3) + (t97 * (-t10 * t46 - t23 * t43) + t93 * (-t10 * t43 + t23 * t46)) * g(2) + (t97 * (-t12 * t46 - t25 * t43) + t93 * (-t12 * t43 + t25 * t46)) * g(1) (-g(1) * t11 - g(2) * t9 - g(3) * t19) * m(7)];
taug  = t1(:);

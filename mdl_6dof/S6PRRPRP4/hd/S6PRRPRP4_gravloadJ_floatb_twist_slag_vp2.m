% Calculate Gravitation load on the joints for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:22
% EndTime: 2019-03-08 21:40:24
% DurationCPUTime: 0.91s
% Computational Cost: add. (443->98), mult. (1107->137), div. (0->0), fcn. (1280->10), ass. (0->56)
t104 = m(7) * pkin(5);
t98 = -m(5) - m(7);
t93 = m(6) - t98;
t97 = -mrSges(6,1) - mrSges(7,1);
t91 = t97 - t104;
t105 = m(7) * (-qJ(6) - pkin(9)) - mrSges(7,3) - mrSges(4,1) + mrSges(5,2);
t102 = -mrSges(4,2) + mrSges(5,3);
t96 = -mrSges(6,2) - mrSges(7,2);
t101 = -m(6) * pkin(9) - pkin(3) * t93 - mrSges(6,3);
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t44 = sin(qJ(5));
t77 = t44 * t45;
t100 = -t102 * t45 - t104 * t77 + t105 * t48 - mrSges(3,1);
t99 = -t101 - t105;
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t67 = cos(pkin(10));
t68 = cos(pkin(6));
t56 = t68 * t67;
t66 = sin(pkin(10));
t26 = t46 * t66 - t49 * t56;
t69 = qJ(4) * t45;
t81 = t26 * t48;
t95 = -pkin(3) * t81 - t26 * t69;
t55 = t68 * t66;
t28 = t46 * t67 + t49 * t55;
t80 = t28 * t48;
t94 = -pkin(3) * t80 - t28 * t69;
t88 = t48 * mrSges(6,3) - t100;
t47 = cos(qJ(5));
t87 = -m(7) * (pkin(5) * t47 + pkin(4)) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t86 = -qJ(4) * t93 + t44 * t91 + t47 * t96 - t102;
t85 = -m(6) * (pkin(4) + pkin(8)) + t87;
t42 = sin(pkin(6));
t79 = t42 * t46;
t78 = t42 * t49;
t76 = t44 * t49;
t75 = t45 * t47;
t74 = t47 * t49;
t70 = pkin(2) * t78 + pkin(8) * t79;
t23 = t26 * pkin(2);
t27 = t46 * t56 + t49 * t66;
t65 = pkin(8) * t27 - t23;
t24 = t28 * pkin(2);
t29 = -t46 * t55 + t49 * t67;
t64 = pkin(8) * t29 - t24;
t63 = t42 * t67;
t62 = t42 * t66;
t31 = t45 * t68 + t48 * t79;
t30 = t45 * t79 - t48 * t68;
t14 = t29 * t48 + t45 * t62;
t13 = t29 * t45 - t48 * t62;
t12 = t27 * t48 - t45 * t63;
t11 = t27 * t45 + t48 * t63;
t1 = [(-m(2) - m(3) - m(4) - t93) * g(3) (-m(4) * t65 - m(6) * (-pkin(9) * t81 - t23 + t95) + t97 * (-t26 * t77 + t27 * t47) + t98 * (t65 + t95) + t96 * (-t26 * t75 - t27 * t44) + t85 * t27 + t88 * t26) * g(2) + (-m(4) * t64 - m(6) * (-pkin(9) * t80 - t24 + t94) + t97 * (-t28 * t77 + t29 * t47) + t96 * (-t28 * t75 - t29 * t44) + t98 * (t64 + t94) + t85 * t29 + t88 * t28) * g(1) + (-m(4) * t70 - t93 * (t69 * t78 + t70) + (t97 * (t45 * t76 + t46 * t47) + t96 * (-t44 * t46 + t45 * t74) + (-m(6) * pkin(4) + t87) * t46 + (t101 * t48 + t100) * t49) * t42) * g(3) (t30 * t99 + t31 * t86) * g(3) + (t11 * t99 + t12 * t86) * g(2) + (t13 * t99 + t14 * t86) * g(1), t93 * (-g(1) * t13 - g(2) * t11 - g(3) * t30) (t96 * (-t30 * t44 + t42 * t74) + t91 * (t30 * t47 + t42 * t76)) * g(3) + (t96 * (-t11 * t44 - t26 * t47) + t91 * (t11 * t47 - t26 * t44)) * g(2) + (t96 * (-t13 * t44 - t28 * t47) + t91 * (t13 * t47 - t28 * t44)) * g(1) (-g(1) * t14 - g(2) * t12 - g(3) * t31) * m(7)];
taug  = t1(:);

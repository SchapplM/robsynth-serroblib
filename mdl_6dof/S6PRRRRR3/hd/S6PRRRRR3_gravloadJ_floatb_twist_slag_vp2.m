% Calculate Gravitation load on the joints for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:58
% EndTime: 2019-03-09 00:47:59
% DurationCPUTime: 1.10s
% Computational Cost: add. (727->110), mult. (1341->155), div. (0->0), fcn. (1576->14), ass. (0->53)
t42 = qJ(4) + qJ(5);
t38 = cos(t42);
t48 = cos(qJ(4));
t40 = t48 * pkin(4);
t27 = pkin(5) * t38 + t40;
t39 = qJ(6) + t42;
t34 = sin(t39);
t35 = cos(t39);
t37 = sin(t42);
t45 = sin(qJ(4));
t96 = -mrSges(4,1) - m(7) * (pkin(3) + t27) - t35 * mrSges(7,1) + t34 * mrSges(7,2) - m(6) * (t40 + pkin(3)) - t38 * mrSges(6,1) + t37 * mrSges(6,2) - m(5) * pkin(3) - t48 * mrSges(5,1) + t45 * mrSges(5,2);
t51 = -pkin(10) - pkin(9);
t95 = mrSges(4,2) + m(7) * (-pkin(11) + t51) - mrSges(7,3) + m(6) * t51 - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t112 = t95 * t46 + t96 * t49 - mrSges(3,1);
t97 = m(4) + m(5) + m(6) + m(7);
t85 = pkin(4) * t45;
t26 = pkin(5) * t37 + t85;
t108 = -m(6) * t85 - m(7) * t26 - t45 * mrSges(5,1) - t37 * mrSges(6,1) - t34 * mrSges(7,1) - t48 * mrSges(5,2) - t38 * mrSges(6,2) - t35 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3);
t107 = pkin(2) * t97 - t112;
t104 = -m(6) * pkin(4) - mrSges(5,1);
t47 = sin(qJ(2));
t78 = cos(pkin(6));
t44 = sin(pkin(6));
t82 = t44 * t49;
t24 = t78 * t46 + t47 * t82;
t50 = cos(qJ(2));
t81 = t44 * t50;
t62 = -t24 * t37 - t38 * t81;
t84 = (-t24 * t34 - t35 * t81) * mrSges(7,1) + (-t24 * t35 + t34 * t81) * mrSges(7,2);
t102 = -t62 * mrSges(6,1) - (-t24 * t38 + t37 * t81) * mrSges(6,2) - t84;
t43 = sin(pkin(12));
t73 = t43 * t78;
t77 = cos(pkin(12));
t22 = -t47 * t73 + t77 * t50;
t16 = t43 * t44 * t46 + t22 * t49;
t21 = t77 * t47 + t50 * t73;
t65 = -t16 * t37 + t21 * t38;
t88 = (-t16 * t34 + t21 * t35) * mrSges(7,1) + (-t16 * t35 - t21 * t34) * mrSges(7,2);
t101 = -t65 * mrSges(6,1) - (-t16 * t38 - t21 * t37) * mrSges(6,2) - t88;
t63 = t78 * t77;
t20 = t43 * t50 + t47 * t63;
t72 = t44 * t77;
t14 = t20 * t49 - t46 * t72;
t19 = t43 * t47 - t50 * t63;
t67 = -t14 * t37 + t19 * t38;
t89 = (-t14 * t34 + t19 * t35) * mrSges(7,1) + (-t14 * t35 - t19 * t34) * mrSges(7,2);
t100 = -t67 * mrSges(6,1) - (-t14 * t38 - t19 * t37) * mrSges(6,2) - t89;
t92 = -t97 * pkin(8) + t108;
t90 = m(7) * pkin(5);
t83 = t44 * t47;
t1 = [(-m(2) - m(3) - t97) * g(3) (-t97 * (pkin(2) * t81 + pkin(8) * t83) + (t108 * t47 + t112 * t50) * t44) * g(3) + (t107 * t19 + t92 * t20) * g(2) + (t107 * t21 + t92 * t22) * g(1) (t95 * t24 + t96 * (-t46 * t83 + t78 * t49)) * g(3) + (t95 * t14 + t96 * (-t20 * t46 - t49 * t72)) * g(2) + (t95 * t16 + t96 * (-t22 * t46 + t43 * t82)) * g(1) (-(-t24 * t48 + t45 * t81) * mrSges(5,2) - m(7) * (-t24 * t26 - t27 * t81) + t104 * (-t24 * t45 - t48 * t81) + t102) * g(3) + (-(-t14 * t48 - t19 * t45) * mrSges(5,2) - m(7) * (-t14 * t26 + t19 * t27) + t104 * (-t14 * t45 + t19 * t48) + t100) * g(2) + (-(-t16 * t48 - t21 * t45) * mrSges(5,2) - m(7) * (-t16 * t26 + t21 * t27) + t104 * (-t16 * t45 + t21 * t48) + t101) * g(1) (-t62 * t90 + t102) * g(3) + (-t67 * t90 + t100) * g(2) + (-t65 * t90 + t101) * g(1), -g(1) * t88 - g(2) * t89 - g(3) * t84];
taug  = t1(:);

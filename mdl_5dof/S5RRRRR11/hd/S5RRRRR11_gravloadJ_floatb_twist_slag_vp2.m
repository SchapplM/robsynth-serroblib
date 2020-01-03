% Calculate Gravitation load on the joints for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:48
% EndTime: 2019-12-31 22:38:52
% DurationCPUTime: 1.02s
% Computational Cost: add. (528->105), mult. (1155->145), div. (0->0), fcn. (1364->12), ass. (0->55)
t46 = cos(qJ(4));
t36 = pkin(4) * t46 + pkin(3);
t40 = qJ(4) + qJ(5);
t37 = sin(t40);
t38 = cos(t40);
t42 = sin(qJ(4));
t91 = m(5) * pkin(3) + m(6) * t36 + t46 * mrSges(5,1) + t38 * mrSges(6,1) - t42 * mrSges(5,2) - t37 * mrSges(6,2) + mrSges(4,1);
t53 = mrSges(4,2) + m(6) * (-pkin(10) - pkin(9)) - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t43 = sin(qJ(3));
t47 = cos(qJ(3));
t105 = t53 * t43 - t91 * t47 - mrSges(3,1);
t104 = m(6) * pkin(4);
t70 = t42 * t104;
t75 = mrSges(3,2) - mrSges(4,3);
t98 = -t37 * mrSges(6,1) - t46 * mrSges(5,2) - t38 * mrSges(6,2);
t103 = -t42 * mrSges(5,1) - t70 + t75 + t98;
t71 = m(4) + m(5) + m(6);
t101 = t71 * pkin(8);
t97 = pkin(2) * t71 - t105;
t94 = mrSges(5,1) + t104;
t90 = -t75 + t101;
t87 = -t101 + t103;
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t48 = cos(qJ(2));
t72 = cos(pkin(5));
t81 = cos(qJ(1));
t62 = t72 * t81;
t24 = t44 * t62 + t45 * t48;
t41 = sin(pkin(5));
t67 = t41 * t81;
t12 = t24 * t47 - t43 * t67;
t23 = t44 * t45 - t48 * t62;
t85 = (-t12 * t37 + t23 * t38) * mrSges(6,1) + (-t12 * t38 - t23 * t37) * mrSges(6,2);
t64 = t45 * t72;
t26 = -t44 * t64 + t81 * t48;
t79 = t41 * t45;
t16 = t26 * t47 + t43 * t79;
t25 = t81 * t44 + t48 * t64;
t5 = -t16 * t37 + t25 * t38;
t6 = t16 * t38 + t25 * t37;
t84 = t5 * mrSges(6,1) - t6 * mrSges(6,2);
t78 = t41 * t47;
t22 = t72 * t43 + t44 * t78;
t77 = t41 * t48;
t82 = (-t22 * t37 - t38 * t77) * mrSges(6,1) + (-t22 * t38 + t37 * t77) * mrSges(6,2);
t80 = t41 * t44;
t73 = t81 * pkin(1) + pkin(7) * t79;
t69 = t26 * pkin(2) + t73;
t65 = -t45 * pkin(1) + pkin(7) * t67;
t11 = -t24 * t43 - t47 * t67;
t7 = -t16 * t42 + t25 * t46;
t15 = t26 * t43 - t45 * t78;
t8 = t16 * t46 + t25 * t42;
t1 = [(-t81 * mrSges(2,1) - m(3) * t73 - t26 * mrSges(3,1) - m(4) * t69 - t16 * mrSges(4,1) - m(5) * (pkin(3) * t16 + t69) - t8 * mrSges(5,1) - t7 * mrSges(5,2) - m(6) * (t16 * t36 + t69) - t6 * mrSges(6,1) - t5 * mrSges(6,2) + (-mrSges(3,3) * t41 + mrSges(2,2)) * t45 + (-t70 - t90) * t25 + t53 * t15) * g(2) + (t45 * mrSges(2,1) + t81 * mrSges(2,2) - m(3) * t65 + t24 * mrSges(3,1) - mrSges(3,3) * t67 + t91 * t12 + (t94 * t42 + t90 - t98) * t23 + t53 * t11 + t71 * (t24 * pkin(2) - t65)) * g(1), (-t71 * (pkin(2) * t77 + pkin(8) * t80) + (t103 * t44 + t105 * t48) * t41) * g(3) + (t97 * t23 + t87 * t24) * g(2) + (t97 * t25 + t87 * t26) * g(1), (t53 * t22 - t91 * (-t43 * t80 + t72 * t47)) * g(3) + (-t91 * t11 + t53 * t12) * g(2) + (t91 * t15 + t53 * t16) * g(1), (-(-t22 * t46 + t42 * t77) * mrSges(5,2) - t82 - t94 * (-t22 * t42 - t46 * t77)) * g(3) + (-(-t12 * t46 - t23 * t42) * mrSges(5,2) - t85 - t94 * (-t12 * t42 + t23 * t46)) * g(2) + (t8 * mrSges(5,2) - t94 * t7 - t84) * g(1), -g(1) * t84 - g(2) * t85 - g(3) * t82];
taug = t1(:);

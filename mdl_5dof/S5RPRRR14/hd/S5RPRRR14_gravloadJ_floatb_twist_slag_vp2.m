% Calculate Gravitation load on the joints for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:54
% EndTime: 2019-12-31 19:16:56
% DurationCPUTime: 0.87s
% Computational Cost: add. (698->98), mult. (1889->150), div. (0->0), fcn. (2380->14), ass. (0->56)
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t105 = m(5) + m(6);
t71 = -pkin(9) * t105 + mrSges(4,2) - mrSges(5,3);
t97 = -t47 * mrSges(6,1) - t50 * mrSges(6,2) + t71;
t52 = cos(qJ(1));
t89 = cos(pkin(11));
t91 = cos(pkin(5));
t75 = t91 * t89;
t87 = sin(pkin(11));
t94 = sin(qJ(1));
t66 = -t52 * t75 + t94 * t87;
t46 = sin(pkin(5));
t88 = sin(pkin(6));
t83 = t46 * t88;
t90 = cos(pkin(6));
t112 = t52 * t83 + t66 * t90;
t73 = t91 * t87;
t34 = t52 * t73 + t89 * t94;
t49 = sin(qJ(3));
t95 = cos(qJ(3));
t18 = t112 * t49 - t34 * t95;
t84 = t46 * t90;
t27 = -t52 * t84 + t66 * t88;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t111 = t18 * t51 - t27 * t48;
t110 = t18 * t48 + t27 * t51;
t103 = m(6) * pkin(4) + t50 * mrSges(6,1) - t47 * mrSges(6,2) + mrSges(5,1);
t82 = -m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t61 = t52 * t87 + t75 * t94;
t108 = t61 * t88 + t94 * t84;
t107 = pkin(3) * t105 + t103 * t51 - t82 * t48 + mrSges(4,1);
t15 = t112 * t95 + t34 * t49;
t101 = t61 * t90 - t94 * t83;
t93 = t46 * t52;
t85 = t46 * t94;
t92 = t52 * pkin(1) + qJ(2) * t85;
t79 = -pkin(1) * t94 + qJ(2) * t93;
t74 = t91 * t88;
t72 = t90 * t89;
t57 = -t34 * pkin(2) - t27 * pkin(8) + t79;
t35 = t52 * t89 - t73 * t94;
t55 = t35 * pkin(2) + t108 * pkin(8) + t92;
t20 = -t101 * t49 + t35 * t95;
t54 = t20 * pkin(3) + t55;
t33 = -t83 * t89 + t90 * t91;
t25 = t49 * t74 + (t49 * t72 + t87 * t95) * t46;
t24 = -t74 * t95 + (t49 * t87 - t72 * t95) * t46;
t19 = t101 * t95 + t35 * t49;
t14 = t25 * t51 + t33 * t48;
t8 = t108 * t48 + t20 * t51;
t7 = -t108 * t51 + t20 * t48;
t2 = t19 * t47 + t50 * t8;
t1 = t19 * t50 - t47 * t8;
t3 = [(-t52 * mrSges(2,1) + t94 * mrSges(2,2) - m(3) * t92 - t35 * mrSges(3,1) + t61 * mrSges(3,2) - mrSges(3,3) * t85 - m(4) * t55 - t20 * mrSges(4,1) - t108 * mrSges(4,3) - m(5) * t54 - t8 * mrSges(5,1) - m(6) * (t8 * pkin(4) + t54) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t82 * t7 + t71 * t19) * g(2) + (t94 * mrSges(2,1) + t52 * mrSges(2,2) - m(3) * t79 + t34 * mrSges(3,1) - t66 * mrSges(3,2) - mrSges(3,3) * t93 - m(4) * t57 - t18 * mrSges(4,1) + t27 * mrSges(4,3) + t82 * t110 - t103 * t111 - t97 * t15 - t105 * (t18 * pkin(3) + t57)) * g(1), (-t91 * g(3) + (-g(1) * t94 + g(2) * t52) * t46) * (t105 + m(3) + m(4)), (t107 * t24 + t97 * t25) * g(3) + (t107 * t15 - t18 * t97) * g(2) + (t107 * t19 + t97 * t20) * g(1), (t82 * t14 - t103 * (-t25 * t48 + t33 * t51)) * g(3) + (-t103 * t110 - t111 * t82) * g(2) + (t103 * t7 + t82 * t8) * g(1), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((t111 * t47 + t15 * t50) * mrSges(6,1) + (t111 * t50 - t15 * t47) * mrSges(6,2)) - g(3) * ((-t14 * t47 + t24 * t50) * mrSges(6,1) + (-t14 * t50 - t24 * t47) * mrSges(6,2))];
taug = t3(:);

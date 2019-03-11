% Calculate Gravitation load on the joints for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:58
% EndTime: 2019-03-09 06:06:59
% DurationCPUTime: 0.69s
% Computational Cost: add. (488->95), mult. (452->104), div. (0->0), fcn. (388->10), ass. (0->54)
t113 = mrSges(6,1) + mrSges(7,1);
t42 = cos(qJ(5));
t114 = t113 * t42;
t112 = mrSges(6,3) + mrSges(7,3);
t104 = mrSges(6,2) + mrSges(7,2);
t35 = pkin(10) + qJ(3);
t32 = qJ(4) + t35;
t26 = sin(t32);
t111 = t26 * t104;
t27 = cos(t32);
t110 = t27 * (-m(6) * pkin(9) - t112);
t108 = t26 * t114;
t107 = t27 * mrSges(5,1) + (-mrSges(5,2) + t112) * t26;
t101 = t27 * pkin(4) + t26 * pkin(9);
t29 = t42 * pkin(5) + pkin(4);
t38 = -qJ(6) - pkin(9);
t103 = -t26 * t38 + t27 * t29;
t106 = -m(6) * t101 - m(7) * t103;
t90 = m(7) * pkin(5);
t52 = -t26 * t29 - t27 * t38;
t88 = pkin(4) * t26;
t30 = sin(t35);
t89 = pkin(3) * t30;
t100 = -m(7) * (t52 - t89) - m(6) * (-t88 - t89) + t108;
t41 = sin(qJ(1));
t43 = cos(qJ(1));
t99 = g(1) * t43 + g(2) * t41;
t98 = m(5) + m(6) + m(7);
t97 = t90 + t113;
t40 = sin(qJ(5));
t96 = -t107 + (t104 * t40 - t114) * t27;
t39 = -pkin(7) - qJ(2);
t95 = -m(3) * qJ(2) + m(4) * t39 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t73 = t41 * t40;
t94 = t41 * t110 - t73 * t111;
t71 = t43 * t40;
t93 = t43 * t110 - t71 * t111;
t92 = m(6) * t88 - m(7) * t52 + t108;
t37 = cos(pkin(10));
t28 = t37 * pkin(2) + pkin(1);
t31 = cos(t35);
t57 = t31 * mrSges(4,1) - t30 * mrSges(4,2);
t91 = m(4) * t28 + mrSges(2,1) + m(3) * pkin(1) + t37 * mrSges(3,1) - sin(pkin(10)) * mrSges(3,2) + t57 + t107;
t25 = pkin(3) * t31;
t72 = t41 * t42;
t70 = t43 * t42;
t54 = mrSges(5,1) * t26 + mrSges(5,2) * t27;
t3 = -t27 * t71 + t72;
t1 = t27 * t73 + t70;
t34 = -pkin(8) + t39;
t10 = t25 + t28;
t4 = t27 * t70 + t73;
t2 = -t27 * t72 + t71;
t5 = [(-t73 * t90 - t98 * (t43 * t10 - t41 * t34) - t113 * t4 - t104 * t3 + t95 * t41 + (t106 - t91) * t43) * g(2) + (-t113 * t2 - t104 * t1 + (t98 * t34 - t90 * t40 + t95) * t43 + (m(5) * t10 - m(6) * (-t10 - t101) - m(7) * (-t10 - t103) + t91) * t41) * g(1) (-g(1) * t41 + g(2) * t43) * (m(3) + m(4) + t98) (t100 * t41 + t94) * g(2) + (t100 * t43 + t93) * g(1) + (-t57 - m(5) * t25 - m(6) * (t25 + t101) - m(7) * (t25 + t103) + t96) * g(3) + (m(5) * t89 + mrSges(4,1) * t30 + mrSges(4,2) * t31 + t54) * t99, t99 * t54 + (t92 * t41 + t94) * g(2) + (t92 * t43 + t93) * g(1) + (t106 + t96) * g(3) (t104 * t42 + t97 * t40) * g(3) * t26 + (t97 * t1 - t104 * t2) * g(2) + (t104 * t4 - t97 * t3) * g(1) (g(3) * t27 - t99 * t26) * m(7)];
taug  = t5(:);

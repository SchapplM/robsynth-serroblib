% Calculate Gravitation load on the joints for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:49
% EndTime: 2019-12-31 21:52:52
% DurationCPUTime: 0.69s
% Computational Cost: add. (314->84), mult. (406->96), div. (0->0), fcn. (356->8), ass. (0->50)
t105 = mrSges(5,1) + mrSges(6,1);
t104 = mrSges(5,3) + mrSges(6,3);
t34 = cos(qJ(4));
t103 = t105 * t34;
t94 = mrSges(5,2) + mrSges(6,2);
t29 = qJ(2) + qJ(3);
t26 = sin(t29);
t102 = t26 * t94;
t27 = cos(t29);
t101 = t27 * (-m(5) * pkin(8) - t104);
t99 = t103 * t26;
t98 = t27 * mrSges(4,1) + (-mrSges(4,2) + t104) * t26;
t92 = t27 * pkin(3) + t26 * pkin(8);
t24 = t34 * pkin(4) + pkin(3);
t30 = -qJ(5) - pkin(8);
t96 = t27 * t24 - t26 * t30;
t97 = -m(5) * t92 - m(6) * t96;
t81 = m(6) * pkin(4);
t45 = -t24 * t26 - t27 * t30;
t79 = pkin(3) * t26;
t32 = sin(qJ(2));
t80 = pkin(2) * t32;
t91 = -m(6) * (t45 - t80) - m(5) * (-t79 - t80) + t99;
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t90 = g(1) * t36 + g(2) * t33;
t89 = m(4) + m(5) + m(6);
t88 = t81 + t105;
t87 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t31 = sin(qJ(4));
t86 = -t98 + (t94 * t31 - t103) * t27;
t66 = t33 * t31;
t85 = t33 * t101 - t66 * t102;
t62 = t36 * t31;
t84 = t36 * t101 - t62 * t102;
t83 = m(5) * t79 - m(6) * t45 + t99;
t35 = cos(qJ(2));
t50 = t35 * mrSges(3,1) - t32 * mrSges(3,2);
t82 = m(3) * pkin(1) + mrSges(2,1) + t50 + t98;
t28 = t35 * pkin(2);
t65 = t33 * t34;
t61 = t36 * t34;
t47 = mrSges(4,1) * t26 + mrSges(4,2) * t27;
t3 = -t27 * t62 + t65;
t1 = t27 * t66 + t61;
t37 = -pkin(7) - pkin(6);
t25 = t28 + pkin(1);
t4 = t27 * t61 + t66;
t2 = -t27 * t65 + t62;
t5 = [(-t66 * t81 - t89 * (t36 * t25 - t33 * t37) - t105 * t4 - t94 * t3 + t87 * t33 + (-t82 + t97) * t36) * g(2) + (-t105 * t2 - t94 * t1 + (-t81 * t31 + t89 * t37 + t87) * t36 + (m(4) * t25 - m(5) * (-t25 - t92) - m(6) * (-t25 - t96) + t82) * t33) * g(1), (t91 * t33 + t85) * g(2) + (t91 * t36 + t84) * g(1) + (-t50 - m(4) * t28 - m(5) * (t28 + t92) - m(6) * (t28 + t96) + t86) * g(3) + (m(4) * t80 + mrSges(3,1) * t32 + mrSges(3,2) * t35 + t47) * t90, t90 * t47 + (t83 * t33 + t85) * g(2) + (t83 * t36 + t84) * g(1) + (t86 + t97) * g(3), (t88 * t31 + t94 * t34) * g(3) * t26 + (t88 * t1 - t94 * t2) * g(2) + (-t88 * t3 + t94 * t4) * g(1), (g(3) * t27 - t90 * t26) * m(6)];
taug = t5(:);

% Calculate Gravitation load on the joints for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:27
% EndTime: 2019-03-09 12:45:29
% DurationCPUTime: 0.94s
% Computational Cost: add. (366->121), mult. (585->132), div. (0->0), fcn. (529->8), ass. (0->54)
t96 = mrSges(6,1) + mrSges(7,1);
t95 = mrSges(6,2) + mrSges(7,2);
t105 = mrSges(3,1) - mrSges(4,2);
t104 = -mrSges(3,2) + mrSges(4,3);
t103 = -mrSges(6,3) - mrSges(7,3);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t92 = g(1) * t43 + g(2) * t40;
t37 = qJ(4) + qJ(5);
t28 = sin(t37);
t38 = sin(qJ(4));
t78 = pkin(4) * t38;
t18 = pkin(5) * t28 + t78;
t29 = cos(t37);
t41 = cos(qJ(4));
t102 = -m(6) * t78 - m(7) * t18 - t38 * mrSges(5,1) - t41 * mrSges(5,2) - t96 * t28 - t95 * t29;
t44 = -pkin(9) - pkin(8);
t36 = -qJ(6) + t44;
t101 = -m(5) * (-pkin(2) - pkin(8)) + mrSges(5,3) - m(6) * (-pkin(2) + t44) - m(7) * (-pkin(2) + t36) - t103;
t97 = -m(6) * pkin(4) - mrSges(5,1);
t42 = cos(qJ(2));
t94 = t95 * t28 * t42;
t39 = sin(qJ(2));
t93 = t104 * t39 + t105 * t42;
t71 = t39 * t40;
t11 = t28 * t43 + t29 * t71;
t12 = -t28 * t71 + t29 * t43;
t91 = -t96 * t11 - t95 * t12;
t70 = t39 * t43;
t10 = t28 * t70 + t29 * t40;
t9 = -t28 * t40 + t29 * t70;
t90 = t95 * t10 - t96 * t9;
t89 = m(4) + m(5) + m(6) + m(7);
t87 = t103 * t42 - t93;
t32 = t41 * pkin(4);
t19 = pkin(5) * t29 + t32;
t85 = -m(5) * pkin(3) - m(6) * (t32 + pkin(3)) - m(7) * (pkin(3) + t19) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t83 = m(7) * pkin(5);
t75 = g(3) * t42;
t33 = t42 * pkin(2);
t72 = t36 * t42;
t69 = t40 * t41;
t68 = t41 * t43;
t65 = t42 * t43;
t64 = t42 * t44;
t30 = t39 * qJ(3);
t63 = t33 + t30;
t62 = t43 * pkin(1) + t40 * pkin(7);
t59 = -pkin(1) - t30;
t13 = -t38 * t40 + t39 * t68;
t15 = t38 * t43 + t39 * t69;
t16 = -t38 * t71 + t68;
t14 = t38 * t70 + t69;
t1 = [(-m(3) * t62 - t14 * mrSges(5,1) - t13 * mrSges(5,2) - t95 * t9 + (-m(5) * pkin(8) - mrSges(5,3)) * t65 - t96 * t10 - t89 * (pkin(2) * t65 + t43 * t30 + t62) + t85 * t40 + (-mrSges(2,1) - m(6) * (t39 * t78 - t64) - m(7) * (t18 * t39 - t72) + t87) * t43) * g(2) + (-t16 * mrSges(5,1) + t15 * mrSges(5,2) - t96 * t12 + t95 * t11 + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t59 - t33) - m(5) * t59 - m(6) * (-pkin(1) + (-qJ(3) - t78) * t39) - m(7) * (-pkin(1) + (-qJ(3) - t18) * t39) + t101 * t42 + t93) * t40 + ((-m(3) - t89) * pkin(7) + t85) * t43) * g(1) (-m(4) * t63 - m(5) * (pkin(8) * t42 + t63) - t42 * mrSges(5,3) - m(6) * (t63 - t64) - m(7) * (t63 - t72) + t102 * t39 + t87) * g(3) + ((m(4) * pkin(2) + t101 + t105) * t39 + (-qJ(3) * t89 + t102 - t104) * t42) * t92 (-t92 * t39 + t75) * t89 -(-mrSges(5,1) * t41 + mrSges(5,2) * t38) * t75 + ((m(6) * t32 + m(7) * t19 + t96 * t29) * t42 - t94) * g(3) + (-t16 * mrSges(5,2) - m(7) * (t18 * t43 + t19 * t71) + t97 * t15 + t91) * g(2) + (t14 * mrSges(5,2) - m(7) * (-t18 * t40 + t19 * t70) + t97 * t13 + t90) * g(1) ((t83 + t96) * t29 * t42 - t94) * g(3) + (-t11 * t83 + t91) * g(2) + (-t9 * t83 + t90) * g(1) (-g(3) * t39 - t92 * t42) * m(7)];
taug  = t1(:);

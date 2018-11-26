% Calculate Gravitation load on the joints for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:07:49
% EndTime: 2018-11-23 17:07:50
% DurationCPUTime: 0.94s
% Computational Cost: add. (373->134), mult. (555->148), div. (0->0), fcn. (498->10), ass. (0->64)
t107 = mrSges(3,1) - mrSges(4,2);
t106 = -mrSges(3,2) + mrSges(4,3);
t105 = -mrSges(6,3) - mrSges(7,3);
t97 = m(6) + m(7);
t104 = t97 + m(4) + m(5);
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t94 = g(1) * t46 + g(2) * t43;
t39 = qJ(4) + pkin(10);
t29 = sin(t39);
t41 = sin(qJ(4));
t84 = pkin(4) * t41;
t18 = pkin(5) * t29 + t84;
t31 = qJ(6) + t39;
t26 = sin(t31);
t27 = cos(t31);
t30 = cos(t39);
t44 = cos(qJ(4));
t103 = -m(6) * t84 - m(7) * t18 - mrSges(5,1) * t41 - t29 * mrSges(6,1) - t26 * mrSges(7,1) - t44 * mrSges(5,2) - t30 * mrSges(6,2) - t27 * mrSges(7,2);
t40 = -qJ(5) - pkin(8);
t38 = -pkin(9) + t40;
t102 = -m(5) * (-pkin(2) - pkin(8)) + mrSges(5,3) - m(6) * (-pkin(2) + t40) - m(7) * (-pkin(2) + t38) - t105;
t96 = -m(6) * pkin(4) - mrSges(5,1);
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t95 = t106 * t42 + t107 * t45;
t91 = t105 * t45 - t95;
t34 = t44 * pkin(4);
t19 = pkin(5) * t30 + t34;
t89 = -m(5) * pkin(3) - m(6) * (t34 + pkin(3)) - m(7) * (pkin(3) + t19) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t76 = t43 * t26;
t77 = t42 * t46;
t5 = t27 * t77 - t76;
t75 = t43 * t27;
t6 = t26 * t77 + t75;
t86 = mrSges(7,1) * t5 - mrSges(7,2) * t6;
t7 = t26 * t46 + t42 * t75;
t8 = t27 * t46 - t42 * t76;
t85 = mrSges(7,1) * t7 + mrSges(7,2) * t8;
t81 = g(3) * t45;
t35 = t45 * pkin(2);
t80 = mrSges(7,1) * t27;
t79 = t38 * t45;
t78 = t40 * t45;
t74 = t43 * t29;
t73 = t43 * t30;
t72 = t43 * t41;
t71 = t43 * t44;
t70 = t44 * t46;
t67 = t45 * t46;
t32 = t42 * qJ(3);
t66 = t35 + t32;
t65 = pkin(1) * t46 + pkin(7) * t43;
t62 = -pkin(1) - t32;
t13 = t42 * t70 - t72;
t15 = t41 * t46 + t42 * t71;
t20 = t45 * t26 * mrSges(7,2);
t16 = -t42 * t72 + t70;
t14 = t41 * t77 + t71;
t12 = t30 * t46 - t42 * t74;
t11 = t29 * t46 + t42 * t73;
t10 = t29 * t77 + t73;
t9 = t30 * t77 - t74;
t1 = [(-m(3) * t65 - t14 * mrSges(5,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t13 * mrSges(5,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (-m(5) * pkin(8) - mrSges(5,3)) * t67 - t104 * (pkin(2) * t67 + t46 * t32 + t65) + t89 * t43 + (-mrSges(2,1) - m(6) * (t42 * t84 - t78) - m(7) * (t18 * t42 - t79) + t91) * t46) * g(2) + (-t16 * mrSges(5,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t15 * mrSges(5,2) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t62 - t35) - m(5) * t62 - m(6) * (-pkin(1) + (-qJ(3) - t84) * t42) - m(7) * (-pkin(1) + (-qJ(3) - t18) * t42) + t102 * t45 + t95) * t43 + ((-m(3) - t104) * pkin(7) + t89) * t46) * g(1) (-m(4) * t66 - m(5) * (pkin(8) * t45 + t66) - t45 * mrSges(5,3) - m(6) * (t66 - t78) - m(7) * (t66 - t79) + t103 * t42 + t91) * g(3) + ((m(4) * pkin(2) + t102 + t107) * t42 + (-qJ(3) * t104 + t103 - t106) * t45) * t94 (-t42 * t94 + t81) * t104, -g(3) * (t20 + (-m(7) * t19 - t80) * t45) + (m(6) * t34 + mrSges(5,1) * t44 + mrSges(6,1) * t30 - mrSges(5,2) * t41 - mrSges(6,2) * t29) * t81 + (-t16 * mrSges(5,2) - t11 * mrSges(6,1) - t12 * mrSges(6,2) - m(7) * (t19 * t42 * t43 + t18 * t46) - t85 + t96 * t15) * g(2) + (t14 * mrSges(5,2) - t9 * mrSges(6,1) + t10 * mrSges(6,2) - m(7) * (-t43 * t18 + t19 * t77) - t86 + t96 * t13) * g(1) (-g(3) * t42 - t45 * t94) * t97, -g(1) * t86 - g(2) * t85 - g(3) * (-t45 * t80 + t20)];
taug  = t1(:);

% Calculate Gravitation load on the joints for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2018-11-23 16:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:29:47
% EndTime: 2018-11-23 16:29:48
% DurationCPUTime: 0.85s
% Computational Cost: add. (369->115), mult. (547->140), div. (0->0), fcn. (513->8), ass. (0->61)
t110 = mrSges(6,1) + mrSges(7,1);
t109 = mrSges(6,2) - mrSges(7,3);
t111 = m(6) + m(7);
t108 = -m(4) - m(5);
t107 = mrSges(6,3) + mrSges(7,2);
t39 = qJ(4) + qJ(5);
t33 = sin(t39);
t34 = cos(t39);
t106 = -t109 * t33 + t110 * t34;
t40 = sin(qJ(4));
t42 = sin(qJ(1));
t43 = cos(qJ(4));
t77 = t42 * t43;
t41 = sin(qJ(3));
t45 = cos(qJ(1));
t78 = t41 * t45;
t17 = t40 * t78 + t77;
t75 = t43 * t45;
t79 = t41 * t42;
t15 = -t40 * t79 + t75;
t95 = -g(1) * t42 + g(2) * t45;
t105 = -m(3) + t108;
t104 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t103 = -m(7) * qJ(6) - mrSges(7,3);
t44 = cos(qJ(3));
t57 = t41 * mrSges(4,1) + t44 * mrSges(4,2);
t102 = mrSges(2,2) - mrSges(3,3) - t57 - m(5) * (pkin(3) * t41 - pkin(8) * t44) + t44 * mrSges(5,3);
t101 = m(7) * pkin(5) + t110;
t100 = mrSges(6,2) + t103;
t32 = pkin(4) * t43 + pkin(3);
t54 = pkin(5) * t34 + qJ(6) * t33;
t98 = -m(7) * (-t32 - t54) + m(6) * t32 + t106;
t46 = -pkin(9) - pkin(8);
t96 = t111 * t46 - t107;
t13 = t33 * t78 + t34 * t42;
t73 = t45 * t34;
t14 = -t33 * t42 + t41 * t73;
t94 = -t109 * t14 - t110 * t13;
t11 = t33 * t79 - t73;
t12 = t33 * t45 + t34 * t79;
t93 = t109 * t12 + t110 * t11;
t91 = -pkin(1) - pkin(7);
t86 = pkin(4) * t40;
t83 = g(3) * t44;
t82 = mrSges(6,2) * t34;
t80 = t40 * t45;
t76 = t42 * t44;
t72 = t45 * t46;
t71 = t17 * pkin(4);
t70 = t45 * pkin(1) + t42 * qJ(2);
t66 = t45 * pkin(7) + t70;
t65 = m(5) * pkin(8) + mrSges(5,3);
t62 = t103 * t34 * t44;
t61 = -t11 * pkin(5) + qJ(6) * t12;
t60 = t13 * pkin(5) - qJ(6) * t14;
t53 = t15 * pkin(4);
t50 = m(5) * pkin(3) + t43 * mrSges(5,1) - t40 * mrSges(5,2);
t36 = t45 * qJ(2);
t18 = -t40 * t42 + t41 * t75;
t16 = t41 * t77 + t80;
t1 = [(-m(3) * t70 - t16 * mrSges(5,1) - t15 * mrSges(5,2) + t107 * t76 + t108 * t66 - t111 * (pkin(4) * t80 + t32 * t79 + t46 * t76 + t66) - t101 * t12 + t100 * t11 + t104 * t45 + t102 * t42) * g(2) + (-t18 * mrSges(5,1) + t17 * mrSges(5,2) - t111 * ((-t86 + t91) * t42 + t32 * t78 + t44 * t72 + t36) - t101 * t14 + t100 * t13 + t105 * t36 + (m(3) * pkin(1) + t108 * t91 - t104) * t42 + (t107 * t44 + t102) * t45) * g(1), t95 * (t111 - t105) ((t107 * t41 + t98 * t44) * t45 - t111 * t41 * t72) * g(2) + (((-m(7) * t54 - t106) * t44 + t96 * t41) * t42 - t111 * t32 * t76) * g(1) + (t57 + (-t65 + t96) * t44 + (t50 + t98) * t41) * g(3) + ((mrSges(4,1) + t50) * t44 + (-mrSges(4,2) + t65) * t41) * t95, -g(3) * ((m(7) * (-pkin(5) * t33 - t86) - t33 * mrSges(7,1)) * t44 - t62) + (m(6) * t86 + mrSges(5,1) * t40 + mrSges(6,1) * t33 + mrSges(5,2) * t43 + t82) * t83 + (-t17 * mrSges(5,1) - t18 * mrSges(5,2) - m(6) * t71 - m(7) * (t60 + t71) + t94) * g(2) + (-t15 * mrSges(5,1) + t16 * mrSges(5,2) - m(6) * t53 - m(7) * (t53 + t61) + t93) * g(1) ((t101 * t33 + t82) * t44 + t62) * g(3) + (-m(7) * t60 + t94) * g(2) + (-m(7) * t61 + t93) * g(1) (-g(1) * t11 + g(2) * t13 - t33 * t83) * m(7)];
taug  = t1(:);

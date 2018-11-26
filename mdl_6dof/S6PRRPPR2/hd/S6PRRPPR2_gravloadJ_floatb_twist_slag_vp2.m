% Calculate Gravitation load on the joints for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:08:39
% EndTime: 2018-11-23 15:08:40
% DurationCPUTime: 0.90s
% Computational Cost: add. (997->109), mult. (1046->148), div. (0->0), fcn. (984->16), ass. (0->59)
t104 = m(6) + m(7);
t110 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t49 = sin(qJ(6));
t52 = cos(qJ(6));
t109 = mrSges(7,1) * t49 + mrSges(7,2) * t52 - mrSges(5,2) + mrSges(6,3);
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t81 = pkin(6) + qJ(2);
t68 = sin(t81) / 0.2e1;
t82 = pkin(6) - qJ(2);
t73 = sin(t82);
t31 = t68 - t73 / 0.2e1;
t45 = sin(pkin(10));
t54 = cos(qJ(2));
t83 = cos(pkin(10));
t63 = t31 * t83 + t45 * t54;
t46 = sin(pkin(6));
t77 = t46 * t83;
t57 = -t50 * t63 - t53 * t77;
t56 = t57 * pkin(3);
t62 = -t31 * t45 + t54 * t83;
t89 = t45 * t46;
t106 = -t50 * t62 + t53 * t89;
t69 = cos(t81) / 0.2e1;
t74 = cos(t82);
t32 = t69 - t74 / 0.2e1;
t47 = cos(pkin(6));
t105 = t32 * t50 + t47 * t53;
t103 = m(5) + t104;
t102 = m(7) * pkin(9) + t110;
t101 = -qJ(5) * t104 - t109;
t44 = qJ(3) + pkin(11);
t42 = sin(t44);
t43 = cos(t44);
t100 = m(4) * pkin(2) + t53 * mrSges(4,1) - t50 * mrSges(4,2) + t109 * t42 + t110 * t43 + mrSges(3,1);
t99 = m(4) * pkin(8) + m(7) * pkin(5) + mrSges(7,1) * t52 - mrSges(7,2) * t49 + mrSges(6,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t51 = sin(qJ(2));
t55 = t74 / 0.2e1 + t69;
t18 = t45 * t51 - t55 * t83;
t98 = t18 * t43;
t21 = t45 * t55 + t51 * t83;
t96 = t21 * t43;
t30 = t68 + t73 / 0.2e1;
t94 = t30 * t43;
t41 = pkin(3) * t53 + pkin(2);
t48 = -qJ(4) - pkin(8);
t87 = -t18 * t41 - t48 * t63;
t86 = -t21 * t41 - t48 * t62;
t85 = t30 * t41 + t32 * t48;
t84 = qJ(5) * t42;
t78 = -pkin(4) * t98 - t18 * t84 + t87;
t76 = -pkin(4) * t96 - t21 * t84 + t86;
t75 = pkin(4) * t94 + t30 * t84 + t85;
t72 = t106 * pkin(3);
t71 = t105 * pkin(3);
t14 = -t32 * t42 - t43 * t47;
t5 = t42 * t62 - t43 * t89;
t3 = t42 * t63 + t43 * t77;
t1 = [(-m(2) - m(3) - m(4) - t103) * g(3) (-m(5) * t85 - m(6) * t75 - m(7) * (pkin(9) * t94 + t75) + t99 * t32 - t100 * t30) * g(3) + (-m(5) * t87 - m(6) * t78 - m(7) * (-pkin(9) * t98 + t78) - t99 * t63 + t100 * t18) * g(2) + (-m(5) * t86 - m(6) * t76 - m(7) * (-pkin(9) * t96 + t76) - t99 * t62 + t100 * t21) * g(1) (-t105 * mrSges(4,1) - (t32 * t53 - t47 * t50) * mrSges(4,2) - m(5) * t71 - t104 * (-pkin(4) * t14 + t71) + t101 * (-t32 * t43 + t42 * t47) + t102 * t14) * g(3) + (-t57 * mrSges(4,1) - (t50 * t77 - t53 * t63) * mrSges(4,2) - m(5) * t56 + t101 * (-t42 * t77 + t43 * t63) + t102 * t3 + t104 * (pkin(4) * t3 - t56)) * g(2) + (-t106 * mrSges(4,1) - (-t50 * t89 - t53 * t62) * mrSges(4,2) - m(5) * t72 - t104 * (-pkin(4) * t5 + t72) + t101 * (t42 * t89 + t43 * t62) + t102 * t5) * g(1), t103 * (-g(1) * t21 - g(2) * t18 + g(3) * t30) t104 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) -g(1) * ((-t21 * t49 + t5 * t52) * mrSges(7,1) + (-t21 * t52 - t49 * t5) * mrSges(7,2)) - g(2) * ((-t18 * t49 + t3 * t52) * mrSges(7,1) + (-t18 * t52 - t3 * t49) * mrSges(7,2)) - g(3) * ((t14 * t52 + t30 * t49) * mrSges(7,1) + (-t14 * t49 + t30 * t52) * mrSges(7,2))];
taug  = t1(:);

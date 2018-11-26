% Calculate Gravitation load on the joints for
% S6PRRRRR2
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:32:38
% EndTime: 2018-11-23 15:32:40
% DurationCPUTime: 1.13s
% Computational Cost: add. (1450->130), mult. (1463->170), div. (0->0), fcn. (1437->18), ass. (0->67)
t69 = qJ(5) + qJ(6);
t65 = sin(t69);
t67 = cos(t69);
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t165 = -t77 * mrSges(6,1) - t67 * mrSges(7,1) + t74 * mrSges(6,2) + t65 * mrSges(7,2) - mrSges(5,1);
t160 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t63 = pkin(5) * t77 + pkin(4);
t70 = qJ(3) + qJ(4);
t66 = sin(t70);
t68 = cos(t70);
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t80 = -pkin(11) - pkin(10);
t163 = -m(4) * pkin(2) - t78 * mrSges(4,1) + t75 * mrSges(4,2) - mrSges(3,1) + (-m(6) * pkin(4) - m(7) * t63 + t165) * t68 + (-m(6) * pkin(10) + m(7) * t80 + t160) * t66;
t150 = -m(7) * pkin(5) - mrSges(6,1);
t71 = sin(pkin(12));
t72 = sin(pkin(6));
t131 = t71 * t72;
t115 = cos(pkin(12));
t114 = pkin(6) - qJ(2);
t103 = sin(t114);
t113 = pkin(6) + qJ(2);
t99 = sin(t113) / 0.2e1;
t55 = t99 - t103 / 0.2e1;
t79 = cos(qJ(2));
t91 = t115 * t79 - t71 * t55;
t158 = t78 * t131 - t91 * t75;
t100 = cos(t113) / 0.2e1;
t104 = cos(t114);
t56 = t100 - t104 / 0.2e1;
t73 = cos(pkin(6));
t157 = t56 * t75 + t73 * t78;
t105 = t72 * t115;
t92 = t115 * t55 + t71 * t79;
t25 = -t68 * t105 - t66 * t92;
t26 = -t66 * t105 + t68 * t92;
t156 = t160 * t26 + t165 * t25;
t27 = t68 * t131 - t66 * t91;
t28 = t66 * t131 + t68 * t91;
t155 = t160 * t28 + t165 * t27;
t40 = t56 * t66 + t68 * t73;
t41 = -t56 * t68 + t66 * t73;
t154 = t160 * t41 + t165 * t40;
t149 = -m(5) - m(6) - m(7);
t143 = -m(4) * pkin(8) - t65 * mrSges(7,1) - t77 * mrSges(6,2) - t67 * mrSges(7,2) + t150 * t74 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t76 = sin(qJ(2));
t84 = t104 / 0.2e1 + t100;
t44 = -t115 * t84 + t71 * t76;
t141 = (-t26 * t65 + t44 * t67) * mrSges(7,1) + (-t26 * t67 - t44 * t65) * mrSges(7,2);
t47 = t115 * t76 + t71 * t84;
t140 = (-t28 * t65 + t47 * t67) * mrSges(7,1) + (-t28 * t67 - t47 * t65) * mrSges(7,2);
t54 = t99 + t103 / 0.2e1;
t139 = (-t41 * t65 - t54 * t67) * mrSges(7,1) + (-t41 * t67 + t54 * t65) * mrSges(7,2);
t137 = t25 * t63 - t26 * t80;
t136 = t27 * t63 - t28 * t80;
t122 = t40 * t63 - t41 * t80;
t108 = t25 * pkin(4) + t26 * pkin(10);
t107 = t27 * pkin(4) + pkin(10) * t28;
t106 = t40 * pkin(4) + pkin(10) * t41;
t102 = t158 * pkin(3);
t101 = t157 * pkin(3);
t87 = -t78 * t105 - t75 * t92;
t85 = t87 * pkin(3);
t81 = -pkin(9) - pkin(8);
t64 = pkin(3) * t78 + pkin(2);
t1 = [(-m(2) - m(3) - m(4) + t149) * g(3) (t149 * (t54 * t64 + t56 * t81) - t143 * t56 + t163 * t54) * g(3) + (t149 * (-t44 * t64 - t81 * t92) + t143 * t92 - t163 * t44) * g(2) + (t149 * (-t47 * t64 - t81 * t91) + t143 * t91 - t163 * t47) * g(1) (-t157 * mrSges(4,1) - (t56 * t78 - t73 * t75) * mrSges(4,2) - m(5) * t101 - m(6) * (t101 + t106) - m(7) * (t101 + t122) + t154) * g(3) + (-t87 * mrSges(4,1) - (t75 * t105 - t78 * t92) * mrSges(4,2) - m(5) * t85 - m(6) * (t108 + t85) - m(7) * (t85 + t137) + t156) * g(2) + (-t158 * mrSges(4,1) - (-t75 * t131 - t78 * t91) * mrSges(4,2) - m(5) * t102 - m(6) * (t102 + t107) - m(7) * (t102 + t136) + t155) * g(1) (-m(6) * t106 - m(7) * t122 + t154) * g(3) + (-m(6) * t108 - m(7) * t137 + t156) * g(2) + (-m(6) * t107 - m(7) * t136 + t155) * g(1) (-(-t41 * t77 + t54 * t74) * mrSges(6,2) - t139 + t150 * (-t41 * t74 - t54 * t77)) * g(3) + (-(-t26 * t77 - t44 * t74) * mrSges(6,2) - t141 + t150 * (-t26 * t74 + t44 * t77)) * g(2) + (-(-t28 * t77 - t47 * t74) * mrSges(6,2) - t140 + t150 * (-t28 * t74 + t47 * t77)) * g(1), -g(1) * t140 - g(2) * t141 - g(3) * t139];
taug  = t1(:);

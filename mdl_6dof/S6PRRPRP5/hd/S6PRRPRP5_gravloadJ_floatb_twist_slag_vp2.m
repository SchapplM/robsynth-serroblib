% Calculate Gravitation load on the joints for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:33
% EndTime: 2019-03-08 21:45:35
% DurationCPUTime: 0.82s
% Computational Cost: add. (475->101), mult. (1202->148), div. (0->0), fcn. (1410->10), ass. (0->56)
t104 = -m(6) - m(7);
t102 = m(5) - t104;
t106 = -mrSges(4,1) + mrSges(5,2);
t105 = mrSges(4,2) - mrSges(5,3);
t103 = -mrSges(6,3) - mrSges(7,2);
t72 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t71 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t101 = t105 * t54 + t106 * t57 - mrSges(3,1);
t100 = -mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t99 = -t103 - t106;
t98 = -t103 * t57 - t101;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t97 = -t102 * qJ(4) - t72 * t53 + t71 * t56 + t105;
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t84 = cos(pkin(10));
t85 = cos(pkin(6));
t66 = t85 * t84;
t83 = sin(pkin(10));
t34 = t55 * t83 - t58 * t66;
t96 = t34 * t57;
t65 = t85 * t83;
t36 = t55 * t84 + t58 * t65;
t95 = t36 * t57;
t52 = sin(pkin(6));
t94 = t52 * t55;
t93 = t52 * t58;
t92 = t53 * t54;
t91 = t53 * t58;
t90 = t54 * t56;
t87 = pkin(2) * t93 + pkin(8) * t94;
t86 = qJ(4) * t54;
t82 = t56 * t93;
t81 = t57 * t93;
t35 = t55 * t66 + t58 * t83;
t80 = -t34 * pkin(2) + pkin(8) * t35;
t37 = -t55 * t65 + t58 * t84;
t79 = -t36 * pkin(2) + t37 * pkin(8);
t75 = t52 * t84;
t74 = t52 * t83;
t73 = pkin(3) * t81 + t86 * t93 + t87;
t68 = -pkin(3) * t96 - t34 * t86 + t80;
t67 = -pkin(3) * t95 - t36 * t86 + t79;
t38 = t54 * t94 - t57 * t85;
t33 = t38 * pkin(3);
t17 = t38 * t56 + t52 * t91;
t15 = t37 * t54 - t57 * t74;
t13 = t35 * t54 + t57 * t75;
t12 = t15 * pkin(3);
t11 = t13 * pkin(3);
t3 = -t15 * t56 + t36 * t53;
t1 = -t13 * t56 + t34 * t53;
t2 = [(-m(2) - m(3) - m(4) - t102) * g(3) (-m(4) * t87 - m(5) * t73 + t103 * t81 + t104 * (pkin(4) * t94 + pkin(9) * t81 + t73) - t71 * (t53 * t94 - t54 * t82) + (-t72 * (t54 * t91 + t55 * t56) + t100 * t55 + t101 * t58) * t52) * g(3) + (-m(4) * t80 - m(5) * t68 + t104 * (t35 * pkin(4) - pkin(9) * t96 + t68) - t72 * (-t34 * t92 + t35 * t56) - t71 * (t34 * t90 + t35 * t53) + t100 * t35 + t98 * t34) * g(2) + (-m(4) * t79 - m(5) * t67 + t104 * (t37 * pkin(4) - pkin(9) * t95 + t67) - t72 * (-t36 * t92 + t37 * t56) - t71 * (t36 * t90 + t37 * t53) + t100 * t37 + t98 * t36) * g(1) (m(5) * t33 + t104 * (-pkin(9) * t38 - t33) + t97 * (t54 * t85 + t57 * t94) + t99 * t38) * g(3) + (m(5) * t11 + t104 * (-pkin(9) * t13 - t11) + t97 * (t35 * t57 - t54 * t75) + t99 * t13) * g(2) + (m(5) * t12 + t104 * (-pkin(9) * t15 - t12) + t97 * (t37 * t57 + t54 * t74) + t99 * t15) * g(1), t102 * (-g(1) * t15 - g(2) * t13 - g(3) * t38) (t71 * (-t38 * t53 + t82) - t72 * t17) * g(3) + (-t71 * (t13 * t53 + t34 * t56) + t72 * t1) * g(2) + (-t71 * (t15 * t53 + t36 * t56) + t72 * t3) * g(1) (-g(1) * t3 - g(2) * t1 + g(3) * t17) * m(7)];
taug  = t2(:);

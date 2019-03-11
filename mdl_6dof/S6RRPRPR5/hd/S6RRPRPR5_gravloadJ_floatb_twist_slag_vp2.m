% Calculate Gravitation load on the joints for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:29:01
% EndTime: 2019-03-09 10:29:04
% DurationCPUTime: 1.30s
% Computational Cost: add. (870->119), mult. (2030->165), div. (0->0), fcn. (2515->14), ass. (0->63)
t47 = pkin(12) + qJ(6);
t45 = sin(t47);
t46 = cos(t47);
t48 = sin(pkin(12));
t50 = cos(pkin(12));
t62 = -m(7) * (pkin(5) * t50 + pkin(4)) - m(6) * pkin(4) - t50 * mrSges(6,1) + t48 * mrSges(6,2) - mrSges(5,1);
t124 = -t46 * mrSges(7,1) + mrSges(7,2) * t45 + t62;
t65 = m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t114 = m(6) + m(7);
t89 = m(5) + t114;
t122 = t50 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t117 = t124 * t55 + t52 * t65 - mrSges(4,1);
t74 = t45 * mrSges(7,1) + t46 * mrSges(7,2);
t115 = m(7) * (pkin(5) * t48 + pkin(9)) + t74 + t48 * mrSges(6,1) + t122 + (m(5) + m(6)) * pkin(9);
t53 = sin(qJ(2));
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t106 = cos(qJ(2));
t92 = cos(pkin(6));
t77 = t92 * t106;
t120 = -t53 * t54 + t56 * t77;
t87 = m(4) + t89;
t90 = sin(pkin(11));
t91 = cos(pkin(11));
t31 = -t106 * t91 + t53 * t90;
t71 = t92 * t90;
t72 = t92 * t91;
t93 = -t106 * t71 - t53 * t72;
t18 = -t31 * t56 + t54 * t93;
t13 = t31 * t54 + t56 * t93;
t108 = pkin(9) * t89 + (m(7) * pkin(5) + mrSges(6,1)) * t48 + t122;
t49 = sin(pkin(6));
t103 = t49 * t54;
t102 = t49 * t56;
t88 = t106 * pkin(2);
t86 = -m(3) * pkin(1) - mrSges(2,1);
t4 = -t102 * t52 - t13 * t55;
t81 = t53 * t92;
t40 = t49 * t88;
t79 = t120 * pkin(2);
t3 = t102 * t55 - t13 * t52;
t29 = -t53 * t56 - t54 * t77;
t32 = -t106 * t90 - t53 * t91;
t63 = t29 * pkin(2);
t60 = mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t49 + t87 * (pkin(2) * t81 + (-pkin(8) - qJ(3)) * t49);
t57 = t106 * t72 - t53 * t71;
t44 = t88 + pkin(1);
t37 = t56 * t44;
t30 = t106 * t56 - t54 * t81;
t28 = -t106 * t54 - t56 * t81;
t25 = t32 * t49;
t24 = t31 * t49;
t20 = -t25 * t55 + t52 * t92;
t19 = -t25 * t52 - t55 * t92;
t17 = t32 * t56 - t54 * t57;
t14 = t32 * t54 + t56 * t57;
t8 = t103 * t52 + t18 * t55;
t7 = -t103 * t55 + t18 * t52;
t2 = -t17 * t45 + t46 * t8;
t1 = -t17 * t46 - t45 * t8;
t5 = [(-t28 * mrSges(3,1) + t120 * mrSges(3,2) + (t44 * t87 - t86) * t54 + t60 * t56 - t65 * t3 - t124 * t4 + (-t108 - t74) * t14 - (pkin(3) * t89 + mrSges(4,1)) * t13) * g(1) + (-m(4) * t37 - t30 * mrSges(3,1) - t18 * mrSges(4,1) - t2 * mrSges(7,1) - t29 * mrSges(3,2) - t1 * mrSges(7,2) + t108 * t17 + t60 * t54 + t86 * t56 + t62 * t8 + t65 * t7 - t89 * (t18 * pkin(3) + t37)) * g(2) (-(mrSges(3,1) * t106 - mrSges(3,2) * t53) * t49 - m(4) * t40 - t89 * (-pkin(3) * t24 + t40) + t115 * t25 - t117 * t24) * g(3) + (-m(4) * t79 - mrSges(3,1) * t120 - mrSges(3,2) * t28 - t89 * (t14 * pkin(3) + t79) + t117 * t14 + t115 * t13) * g(2) + (-m(4) * t63 - mrSges(3,1) * t29 + mrSges(3,2) * t30 - t89 * (pkin(3) * t17 + t63) + t117 * t17 - t115 * t18) * g(1) ((-g(1) * t54 + g(2) * t56) * t49 - g(3) * t92) * t87 (-t124 * t19 + t65 * t20) * g(3) + (-t124 * t3 + t65 * t4) * g(2) + (-t124 * t7 + t65 * t8) * g(1), t114 * (-g(1) * t7 - g(2) * t3 - g(3) * t19) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t46 - t4 * t45) * mrSges(7,1) + (t14 * t45 - t4 * t46) * mrSges(7,2)) - g(3) * ((-t20 * t45 + t24 * t46) * mrSges(7,1) + (-t20 * t46 - t24 * t45) * mrSges(7,2))];
taug  = t5(:);

% Calculate Gravitation load on the joints for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:11
% EndTime: 2019-03-08 21:18:14
% DurationCPUTime: 0.93s
% Computational Cost: add. (442->90), mult. (1026->119), div. (0->0), fcn. (1180->12), ass. (0->51)
t32 = pkin(11) + qJ(6);
t30 = sin(t32);
t31 = cos(t32);
t33 = sin(pkin(11));
t35 = cos(pkin(11));
t104 = -t30 * mrSges(7,1) - t35 * mrSges(6,2) - t31 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t33;
t102 = m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2);
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t103 = t102 * t39 + t104 * t37 - mrSges(3,1);
t89 = m(6) + m(7);
t86 = m(5) + t89;
t66 = qJ(4) * t37;
t99 = -pkin(3) * t39 - t66;
t91 = pkin(3) * t86 + mrSges(6,3) - t102;
t90 = -t35 * mrSges(6,1) - t31 * mrSges(7,1) + t33 * mrSges(6,2) + t30 * mrSges(7,2) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t38 = sin(qJ(2));
t40 = cos(qJ(2));
t63 = cos(pkin(10));
t64 = cos(pkin(6));
t49 = t64 * t63;
t62 = sin(pkin(10));
t14 = t38 * t62 - t40 * t49;
t88 = t99 * t14;
t48 = t64 * t62;
t16 = t38 * t63 + t40 * t48;
t87 = t99 * t16;
t82 = -t86 * qJ(4) + t104;
t29 = pkin(5) * t35 + pkin(4);
t79 = -m(6) * (pkin(4) + pkin(8)) - m(7) * (pkin(8) + t29) + t90;
t78 = t39 * mrSges(6,3) - t103;
t34 = sin(pkin(6));
t73 = t34 * t38;
t72 = t34 * t40;
t68 = t39 * t40;
t67 = pkin(2) * t72 + pkin(8) * t73;
t11 = t14 * pkin(2);
t15 = t38 * t49 + t40 * t62;
t58 = pkin(8) * t15 - t11;
t12 = t16 * pkin(2);
t17 = -t38 * t48 + t40 * t63;
t57 = pkin(8) * t17 - t12;
t56 = t34 * t63;
t55 = t34 * t62;
t19 = t37 * t64 + t39 * t73;
t18 = t37 * t73 - t39 * t64;
t6 = t17 * t39 + t37 * t55;
t5 = t17 * t37 - t39 * t55;
t4 = t15 * t39 - t37 * t56;
t3 = t15 * t37 + t39 * t56;
t1 = [(-m(2) - m(3) - m(4) - t86) * g(3) (-m(4) * t67 - t86 * (t34 * pkin(3) * t68 + t66 * t72 + t67) + (-mrSges(6,3) * t68 + (-m(6) * pkin(4) - m(7) * t29 + t90) * t38 + t103 * t40) * t34) * g(3) + (-m(4) * t58 - m(5) * (t58 + t88) - t89 * (-t11 + t88) + t79 * t15 + t78 * t14) * g(2) + (-m(4) * t57 - m(5) * (t57 + t87) - t89 * (-t12 + t87) + t79 * t17 + t78 * t16) * g(1) (t91 * t18 + t82 * t19) * g(3) + (t91 * t3 + t82 * t4) * g(2) + (t91 * t5 + t82 * t6) * g(1), t86 * (-g(1) * t5 - g(2) * t3 - g(3) * t18) t89 * (-g(1) * t6 - g(2) * t4 - g(3) * t19) -g(1) * ((-t16 * t30 + t31 * t5) * mrSges(7,1) + (-t16 * t31 - t30 * t5) * mrSges(7,2)) - g(2) * ((-t14 * t30 + t3 * t31) * mrSges(7,1) + (-t14 * t31 - t3 * t30) * mrSges(7,2)) - g(3) * ((t18 * t31 + t30 * t72) * mrSges(7,1) + (-t18 * t30 + t31 * t72) * mrSges(7,2))];
taug  = t1(:);

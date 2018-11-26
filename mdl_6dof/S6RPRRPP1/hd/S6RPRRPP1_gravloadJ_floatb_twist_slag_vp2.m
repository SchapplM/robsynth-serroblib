% Calculate Gravitation load on the joints for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2018-11-23 16:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:11:16
% EndTime: 2018-11-23 16:11:16
% DurationCPUTime: 0.74s
% Computational Cost: add. (474->95), mult. (477->113), div. (0->0), fcn. (439->10), ass. (0->53)
t87 = mrSges(6,1) + mrSges(7,1);
t86 = -mrSges(6,2) + mrSges(7,3);
t72 = m(6) + m(7);
t81 = pkin(4) * t72 + mrSges(5,1);
t85 = mrSges(4,2) - mrSges(7,2) - mrSges(6,3);
t26 = qJ(4) + pkin(10);
t21 = sin(t26);
t23 = cos(t26);
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t84 = m(5) * pkin(3) + t32 * mrSges(5,1) - t29 * mrSges(5,2) + t86 * t21 + t23 * t87;
t27 = qJ(1) + pkin(9);
t22 = sin(t27);
t24 = cos(t27);
t33 = cos(qJ(3));
t57 = t33 * t29;
t7 = t22 * t32 - t24 * t57;
t82 = g(1) * t24 + g(2) * t22;
t20 = pkin(4) * t32 + pkin(3);
t16 = t33 * t20;
t28 = -qJ(5) - pkin(8);
t30 = sin(qJ(3));
t63 = t28 * t30;
t80 = t72 * (t16 - t63);
t79 = -m(4) - m(5);
t78 = mrSges(3,2) - mrSges(4,3);
t77 = -t33 * mrSges(4,1) + t85 * t30;
t75 = m(7) * pkin(5) + t87;
t74 = m(7) * qJ(6) + t86;
t73 = t30 * mrSges(5,3) + mrSges(3,1) - t77;
t31 = sin(qJ(1));
t71 = pkin(1) * t31;
t67 = g(3) * t30;
t34 = cos(qJ(1));
t25 = t34 * pkin(1);
t66 = t22 * t29;
t64 = t24 * t29;
t59 = t33 * t21;
t58 = t33 * t23;
t56 = t33 * t32;
t54 = t24 * pkin(2) + t22 * pkin(7) + t25;
t53 = m(5) * pkin(8) + mrSges(5,3);
t52 = t24 * pkin(7) - t71;
t49 = pkin(3) * t33 + pkin(8) * t30;
t43 = pkin(5) * t23 + qJ(6) * t21;
t5 = t22 * t57 + t24 * t32;
t8 = t24 * t56 + t66;
t6 = -t22 * t56 + t64;
t4 = t21 * t22 + t24 * t58;
t3 = -t22 * t23 + t24 * t59;
t2 = -t24 * t21 + t22 * t58;
t1 = t22 * t59 + t23 * t24;
t9 = [(-m(3) * t25 - t34 * mrSges(2,1) - t8 * mrSges(5,1) + t31 * mrSges(2,2) - t7 * mrSges(5,2) + t79 * t54 - t72 * (pkin(4) * t66 + t54) - t75 * t4 - t74 * t3 + t78 * t22 + (-m(5) * t49 - t73 - t80) * t24) * g(2) + (m(3) * t71 + t31 * mrSges(2,1) - t6 * mrSges(5,1) + t34 * mrSges(2,2) - t5 * mrSges(5,2) + t79 * t52 - t72 * (pkin(4) * t64 + t22 * t63 + t52) + t78 * t24 + t75 * t2 + t74 * t1 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t49) - t72 * (-pkin(2) - t16) + t73) * t22) * g(1) (-m(3) - t72 + t79) * g(3) (-t80 + t77) * g(3) + ((-m(7) * t43 - t84) * g(3) + t82 * (t72 * t28 - t53 + t85)) * t33 + (-t53 * g(3) + t82 * (mrSges(4,1) + m(6) * t20 - m(7) * (-t20 - t43) + t84)) * t30 (mrSges(5,2) * t32 + t75 * t21 - t74 * t23 + t81 * t29) * t67 + (-t6 * mrSges(5,2) + t75 * t1 - t74 * t2 + t81 * t5) * g(2) + (t8 * mrSges(5,2) + t75 * t3 - t74 * t4 - t81 * t7) * g(1) (t33 * g(3) - t30 * t82) * t72 (-g(1) * t3 - g(2) * t1 - t21 * t67) * m(7)];
taug  = t9(:);

% Calculate Gravitation load on the joints for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2018-11-23 16:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:13:36
% EndTime: 2018-11-23 16:13:37
% DurationCPUTime: 0.76s
% Computational Cost: add. (414->102), mult. (550->113), div. (0->0), fcn. (521->8), ass. (0->52)
t87 = mrSges(5,1) + mrSges(6,1);
t76 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t86 = -mrSges(5,3) - mrSges(6,2);
t85 = mrSges(7,3) + t86;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t84 = -t76 * t28 + t87 * t30;
t29 = sin(qJ(1));
t83 = g(2) * t29;
t24 = pkin(9) + qJ(3);
t22 = sin(t24);
t67 = g(3) * t22;
t26 = cos(pkin(9));
t82 = mrSges(2,1) + m(3) * pkin(1) + t26 * mrSges(3,1) - sin(pkin(9)) * mrSges(3,2);
t81 = m(6) + m(7);
t31 = cos(qJ(1));
t79 = g(1) * t31 + t83;
t78 = -m(5) - t81;
t77 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t23 = cos(t24);
t45 = m(7) * (-pkin(4) - pkin(5)) - mrSges(7,1);
t55 = qJ(5) * t28;
t49 = -pkin(3) - t55;
t70 = pkin(4) * t30;
t75 = (m(7) * qJ(6) + t85) * t23 + (-m(7) * t49 - t30 * t45 - m(6) * (t49 - t70) + m(5) * pkin(3) + t84) * t22;
t44 = t23 * mrSges(4,1) - t22 * mrSges(4,2);
t74 = t86 * t22 - t44;
t73 = m(7) * pkin(5) + mrSges(7,1);
t72 = pkin(8) * t78;
t71 = t73 + t87;
t17 = t22 * pkin(8);
t18 = t23 * pkin(3);
t64 = t22 * t31;
t63 = t23 * t31;
t27 = -pkin(7) - qJ(2);
t62 = t27 * t31;
t60 = t29 * t28;
t59 = t29 * t30;
t58 = t30 * t31;
t57 = t31 * t28;
t56 = t18 + t17;
t54 = qJ(6) * t22;
t19 = pkin(2) * t26 + pkin(1);
t53 = -t19 - t18;
t48 = t31 * t19 - t29 * t27;
t47 = t56 + (t55 + t70) * t23;
t40 = pkin(3) * t63 + pkin(8) * t64 + t48;
t8 = t23 * t58 + t60;
t7 = t23 * t57 - t59;
t6 = t23 * t59 - t57;
t5 = t23 * t60 + t58;
t1 = [(-m(4) * t48 - m(5) * t40 - t81 * (t8 * pkin(4) + t7 * qJ(5) + t40) - t71 * t8 + t76 * t7 + t85 * t64 + (m(7) * t54 - t44 - t82) * t31 + t77 * t29) * g(2) + (m(5) * t62 - t81 * (-t6 * pkin(4) - t5 * qJ(5) - t62) + t71 * t6 - t76 * t5 + (m(4) * t27 + t77) * t31 + (m(4) * t19 - m(7) * t53 - (m(7) * (-pkin(8) + qJ(6)) + mrSges(7,3)) * t22 + (-m(5) - m(6)) * (t53 - t17) - t74 + t82) * t29) * g(1) (-g(1) * t29 + g(2) * t31) * (m(3) + m(4) - t78) t79 * (mrSges(4,1) * t22 + mrSges(4,2) * t23) + (t75 * t31 + t63 * t72) * g(1) + (-m(5) * t56 - m(6) * t47 - m(7) * (t47 - t54) + t22 * mrSges(7,3) + (-t73 * t30 - t84) * t23 + t74) * g(3) + (t23 * t72 + t75) * t83 (-t81 * (-t5 * pkin(4) + qJ(5) * t6) + t76 * t6 + t71 * t5) * g(2) + (-t81 * (-t7 * pkin(4) + qJ(5) * t8) + t76 * t8 + t71 * t7) * g(1) + ((-t81 * qJ(5) + t76) * t30 + (m(6) * pkin(4) - t45 + t87) * t28) * t67, t81 * (-g(1) * t7 - g(2) * t5 - t28 * t67) (-g(3) * t23 + t79 * t22) * m(7)];
taug  = t1(:);

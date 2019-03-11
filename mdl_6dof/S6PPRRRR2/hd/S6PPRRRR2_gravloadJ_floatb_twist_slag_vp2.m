% Calculate Gravitation load on the joints for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:16
% EndTime: 2019-03-08 19:03:17
% DurationCPUTime: 0.74s
% Computational Cost: add. (897->85), mult. (2354->132), div. (0->0), fcn. (2987->16), ass. (0->55)
t90 = m(5) + m(6) + m(7);
t36 = qJ(5) + qJ(6);
t34 = sin(t36);
t35 = cos(t36);
t37 = sin(qJ(5));
t40 = cos(qJ(5));
t96 = mrSges(5,1) + m(7) * (pkin(5) * t40 + pkin(4)) + t35 * mrSges(7,1) - t34 * mrSges(7,2) + m(6) * pkin(4) + t40 * mrSges(6,1) - t37 * mrSges(6,2);
t88 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) - mrSges(6,3);
t94 = -m(7) * pkin(5) - mrSges(6,1);
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t95 = pkin(3) * t90 - t88 * t38 + t96 * t41 + mrSges(4,1);
t71 = sin(pkin(13));
t72 = sin(pkin(12));
t56 = t72 * t71;
t75 = cos(pkin(13));
t76 = cos(pkin(12));
t63 = t76 * t75;
t78 = cos(pkin(6));
t47 = -t78 * t63 + t56;
t73 = sin(pkin(7));
t74 = sin(pkin(6));
t60 = t74 * t73;
t77 = cos(pkin(7));
t93 = t47 * t77 + t76 * t60;
t57 = t72 * t75;
t61 = t76 * t71;
t48 = t78 * t57 + t61;
t59 = t74 * t72;
t92 = t48 * t77 - t73 * t59;
t91 = t75 * t77 * t74 + t78 * t73;
t85 = m(3) + m(4) + t90;
t84 = -t34 * mrSges(7,1) - t40 * mrSges(6,2) - t35 * mrSges(7,2) - t90 * pkin(9) + t94 * t37 + mrSges(4,2) - mrSges(5,3);
t27 = t78 * t61 + t57;
t39 = sin(qJ(3));
t79 = cos(qJ(3));
t13 = t27 * t39 + t93 * t79;
t14 = t27 * t79 - t93 * t39;
t62 = t76 * t74;
t22 = t47 * t73 - t77 * t62;
t8 = t14 * t41 + t22 * t38;
t82 = (t13 * t35 - t34 * t8) * mrSges(7,1) + (-t13 * t34 - t35 * t8) * mrSges(7,2);
t28 = -t78 * t56 + t63;
t16 = t28 * t79 - t92 * t39;
t23 = t48 * t73 + t77 * t59;
t10 = t16 * t41 + t23 * t38;
t15 = t28 * t39 + t92 * t79;
t81 = (-t10 * t34 + t15 * t35) * mrSges(7,1) + (-t10 * t35 - t15 * t34) * mrSges(7,2);
t58 = t74 * t71;
t21 = t91 * t39 + t79 * t58;
t26 = -t75 * t60 + t78 * t77;
t18 = t21 * t41 + t26 * t38;
t20 = t39 * t58 - t91 * t79;
t80 = (-t18 * t34 + t20 * t35) * mrSges(7,1) + (-t18 * t35 - t20 * t34) * mrSges(7,2);
t1 = [(-m(2) - t85) * g(3) (-t59 * g(1) + t62 * g(2) - t78 * g(3)) * t85 (t95 * t20 + t84 * t21) * g(3) + (t95 * t13 + t84 * t14) * g(2) + (t95 * t15 + t84 * t16) * g(1) (t88 * t18 - t96 * (-t21 * t38 + t26 * t41)) * g(3) + (t88 * t8 - t96 * (-t14 * t38 + t22 * t41)) * g(2) + (-t96 * (-t16 * t38 + t23 * t41) + t88 * t10) * g(1) (-(-t18 * t40 - t20 * t37) * mrSges(6,2) - t80 + t94 * (-t18 * t37 + t20 * t40)) * g(3) + (-(-t13 * t37 - t40 * t8) * mrSges(6,2) - t82 + t94 * (t13 * t40 - t37 * t8)) * g(2) + (-(-t10 * t40 - t15 * t37) * mrSges(6,2) - t81 + t94 * (-t10 * t37 + t15 * t40)) * g(1), -g(1) * t81 - g(2) * t82 - g(3) * t80];
taug  = t1(:);

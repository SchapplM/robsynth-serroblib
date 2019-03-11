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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:18
% EndTime: 2019-03-08 21:04:21
% DurationCPUTime: 1.00s
% Computational Cost: add. (493->107), mult. (878->153), div. (0->0), fcn. (984->12), ass. (0->57)
t98 = m(6) + m(7);
t106 = -mrSges(5,1) + mrSges(6,2);
t105 = mrSges(5,2) - mrSges(6,3);
t39 = sin(qJ(6));
t42 = cos(qJ(6));
t88 = -m(4) * pkin(8) - m(7) * pkin(5) - t42 * mrSges(7,1) + t39 * mrSges(7,2) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t35 = qJ(3) + pkin(11);
t33 = sin(t35);
t34 = cos(t35);
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t104 = m(4) * pkin(2) + t43 * mrSges(4,1) - t40 * mrSges(4,2) - t105 * t33 - t106 * t34 + mrSges(3,1);
t103 = t39 * mrSges(7,1) + t42 * mrSges(7,2);
t36 = sin(pkin(10));
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t68 = cos(pkin(10));
t69 = cos(pkin(6));
t52 = t69 * t68;
t19 = t36 * t44 + t41 * t52;
t37 = sin(pkin(6));
t61 = t37 * t68;
t47 = -t19 * t40 - t43 * t61;
t45 = t47 * pkin(3);
t77 = t37 * t41;
t100 = -t40 * t77 + t69 * t43;
t62 = t36 * t69;
t21 = -t41 * t62 + t68 * t44;
t76 = t37 * t43;
t99 = -t21 * t40 + t36 * t76;
t94 = m(5) + t98;
t93 = -m(7) * pkin(9) - mrSges(7,3);
t92 = -t93 - t106;
t91 = -t98 * qJ(5) - t103 + t105;
t89 = t34 * mrSges(7,3) + t103 * t33 + t104;
t18 = t36 * t41 - t44 * t52;
t86 = t18 * t34;
t20 = t68 * t41 + t44 * t62;
t84 = t20 * t34;
t78 = t36 * t37;
t75 = t37 * t44;
t74 = t39 * t44;
t73 = t42 * t44;
t32 = pkin(3) * t43 + pkin(2);
t38 = -qJ(4) - pkin(8);
t72 = -t18 * t32 - t19 * t38;
t71 = -t20 * t32 - t21 * t38;
t70 = qJ(5) * t33;
t63 = -pkin(4) * t86 - t18 * t70 + t72;
t59 = -pkin(4) * t84 - t20 * t70 + t71;
t58 = t99 * pkin(3);
t53 = t100 * pkin(3);
t24 = t32 * t75;
t14 = t33 * t77 - t69 * t34;
t5 = t21 * t33 - t34 * t78;
t3 = t19 * t33 + t34 * t61;
t1 = [(-m(2) - m(3) - m(4) - t94) * g(3) (-m(5) * t72 - m(6) * t63 - m(7) * (-pkin(9) * t86 + t63) + t88 * t19 + t89 * t18) * g(2) + (-m(5) * t71 - m(6) * t59 - m(7) * (-pkin(9) * t84 + t59) + t88 * t21 + t89 * t20) * g(1) + (-m(5) * t24 - t98 * (t24 + (pkin(4) * t34 + t70) * t75) + ((-t74 * mrSges(7,1) - t73 * mrSges(7,2)) * t33 + (t94 * t38 + t88) * t41 + (t93 * t34 - t104) * t44) * t37) * g(3) (-t100 * mrSges(4,1) - (-t69 * t40 - t41 * t76) * mrSges(4,2) - m(5) * t53 - t98 * (-t14 * pkin(4) + t53) + t91 * (t69 * t33 + t34 * t77) + t92 * t14) * g(3) + (-t47 * mrSges(4,1) - (-t19 * t43 + t40 * t61) * mrSges(4,2) - m(5) * t45 + t91 * (t19 * t34 - t33 * t61) + t92 * t3 + t98 * (t3 * pkin(4) - t45)) * g(2) + (-t99 * mrSges(4,1) - (-t21 * t43 - t40 * t78) * mrSges(4,2) - m(5) * t58 - t98 * (-t5 * pkin(4) + t58) + t91 * (t21 * t34 + t33 * t78) + t92 * t5) * g(1), t94 * (-g(1) * t20 - g(2) * t18 + g(3) * t75) t98 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) -g(1) * ((-t20 * t39 + t42 * t5) * mrSges(7,1) + (-t20 * t42 - t39 * t5) * mrSges(7,2)) - g(2) * ((-t18 * t39 + t3 * t42) * mrSges(7,1) + (-t18 * t42 - t3 * t39) * mrSges(7,2)) - g(3) * ((t14 * t42 + t37 * t74) * mrSges(7,1) + (-t14 * t39 + t37 * t73) * mrSges(7,2))];
taug  = t1(:);

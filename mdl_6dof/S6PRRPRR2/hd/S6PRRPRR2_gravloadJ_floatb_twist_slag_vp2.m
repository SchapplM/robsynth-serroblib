% Calculate Gravitation load on the joints for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:55:59
% EndTime: 2019-03-08 21:56:03
% DurationCPUTime: 1.20s
% Computational Cost: add. (619->101), mult. (1034->146), div. (0->0), fcn. (1182->14), ass. (0->45)
t34 = qJ(5) + qJ(6);
t31 = sin(t34);
t32 = cos(t34);
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t108 = mrSges(5,1) + m(7) * (pkin(5) * t41 + pkin(4)) + t32 * mrSges(7,1) - t31 * mrSges(7,2) + m(6) * pkin(4) + t41 * mrSges(6,1) - t38 * mrSges(6,2);
t88 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) - mrSges(6,3);
t33 = qJ(3) + pkin(12);
t29 = sin(t33);
t30 = cos(t33);
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t90 = m(4) * pkin(2) + t42 * mrSges(4,1) - t39 * mrSges(4,2) + t108 * t30 - t88 * t29 + mrSges(3,1);
t98 = -m(7) * pkin(5) - mrSges(6,1);
t92 = m(5) + m(6) + m(7);
t73 = cos(pkin(6));
t36 = sin(pkin(6));
t40 = sin(qJ(2));
t78 = t36 * t40;
t103 = -t39 * t78 + t73 * t42;
t43 = cos(qJ(2));
t35 = sin(pkin(11));
t66 = t35 * t73;
t72 = cos(pkin(11));
t20 = -t40 * t66 + t72 * t43;
t77 = t36 * t42;
t102 = -t20 * t39 + t35 * t77;
t89 = -m(4) * pkin(8) - t31 * mrSges(7,1) - t41 * mrSges(6,2) - t32 * mrSges(7,2) + t98 * t38 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t55 = t73 * t72;
t18 = t35 * t43 + t40 * t55;
t65 = t36 * t72;
t94 = -t18 * t39 - t42 * t65;
t17 = t35 * t40 - t43 * t55;
t8 = t18 * t30 - t29 * t65;
t86 = (t17 * t32 - t31 * t8) * mrSges(7,1) + (-t17 * t31 - t32 * t8) * mrSges(7,2);
t79 = t35 * t36;
t10 = t20 * t30 + t29 * t79;
t19 = t72 * t40 + t43 * t66;
t85 = (-t10 * t31 + t19 * t32) * mrSges(7,1) + (-t10 * t32 - t19 * t31) * mrSges(7,2);
t14 = t73 * t29 + t30 * t78;
t76 = t36 * t43;
t84 = (-t14 * t31 - t32 * t76) * mrSges(7,1) + (-t14 * t32 + t31 * t76) * mrSges(7,2);
t37 = -qJ(4) - pkin(8);
t28 = pkin(3) * t42 + pkin(2);
t1 = [(-m(2) - m(3) - m(4) - t92) * g(3) (-t92 * (-t17 * t28 - t18 * t37) + t89 * t18 + t90 * t17) * g(2) + (-t92 * (-t19 * t28 - t20 * t37) + t89 * t20 + t90 * t19) * g(1) + (-t92 * t28 * t76 + (-t90 * t43 + (t92 * t37 + t89) * t40) * t36) * g(3) (-t103 * mrSges(4,1) - (-t73 * t39 - t40 * t77) * mrSges(4,2) + t88 * t14 - t108 * (-t29 * t78 + t73 * t30)) * g(3) + (-(-t18 * t42 + t39 * t65) * mrSges(4,2) - mrSges(4,1) * t94 + t88 * t8 - t108 * (-t18 * t29 - t30 * t65)) * g(2) + (-t102 * mrSges(4,1) - (-t20 * t42 - t39 * t79) * mrSges(4,2) - t108 * (-t20 * t29 + t30 * t79) + t88 * t10) * g(1) + (-g(1) * t102 - t94 * g(2) - g(3) * t103) * t92 * pkin(3), t92 * (-g(1) * t19 - g(2) * t17 + g(3) * t76) (-(-t14 * t41 + t38 * t76) * mrSges(6,2) - t84 + t98 * (-t14 * t38 - t41 * t76)) * g(3) + (-(-t17 * t38 - t41 * t8) * mrSges(6,2) - t86 + t98 * (t17 * t41 - t38 * t8)) * g(2) + (-(-t10 * t41 - t19 * t38) * mrSges(6,2) - t85 + t98 * (-t10 * t38 + t19 * t41)) * g(1), -g(1) * t85 - g(2) * t86 - g(3) * t84];
taug  = t1(:);

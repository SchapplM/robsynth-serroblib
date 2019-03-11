% Calculate Gravitation load on the joints for
% S6PRRPRR5
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:24
% EndTime: 2019-03-08 22:16:26
% DurationCPUTime: 0.90s
% Computational Cost: add. (599->89), mult. (1129->124), div. (0->0), fcn. (1311->14), ass. (0->48)
t36 = pkin(12) + qJ(5);
t33 = qJ(6) + t36;
t28 = sin(t33);
t29 = cos(t33);
t39 = cos(pkin(12));
t30 = t39 * pkin(4) + pkin(3);
t31 = sin(t36);
t32 = cos(t36);
t37 = sin(pkin(12));
t84 = mrSges(4,1) + m(7) * (pkin(5) * t32 + t30) + t29 * mrSges(7,1) - t28 * mrSges(7,2) + m(6) * t30 + t32 * mrSges(6,1) - t31 * mrSges(6,2) + m(5) * pkin(3) + t39 * mrSges(5,1) - t37 * mrSges(5,2);
t40 = -pkin(9) - qJ(4);
t83 = mrSges(4,2) - m(5) * qJ(4) - mrSges(5,3) + m(7) * (-pkin(10) + t40) - mrSges(7,3) + m(6) * t40 - mrSges(6,3);
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t97 = t83 * t41 - t84 * t43 - mrSges(3,1);
t88 = m(5) + m(6) + m(7);
t85 = m(4) + t88;
t75 = pkin(4) * t37;
t94 = -m(6) * t75 - m(7) * (pkin(5) * t31 + t75) - t37 * mrSges(5,1) - t31 * mrSges(6,1) - t28 * mrSges(7,1) - t39 * mrSges(5,2) - t32 * mrSges(6,2) - t29 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3);
t93 = pkin(2) * t85 - t97;
t90 = -m(7) * pkin(5) - mrSges(6,1);
t80 = -t85 * pkin(8) + t94;
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t70 = cos(pkin(11));
t71 = cos(pkin(6));
t57 = t71 * t70;
t69 = sin(pkin(11));
t13 = t69 * t42 - t44 * t57;
t14 = t42 * t57 + t69 * t44;
t38 = sin(pkin(6));
t65 = t38 * t70;
t8 = t14 * t43 - t41 * t65;
t78 = (t13 * t29 - t28 * t8) * mrSges(7,1) + (-t13 * t28 - t29 * t8) * mrSges(7,2);
t56 = t71 * t69;
t16 = -t42 * t56 + t70 * t44;
t64 = t38 * t69;
t10 = t16 * t43 + t41 * t64;
t15 = t70 * t42 + t44 * t56;
t77 = (-t10 * t28 + t15 * t29) * mrSges(7,1) + (-t10 * t29 - t15 * t28) * mrSges(7,2);
t74 = t38 * t42;
t18 = t71 * t41 + t43 * t74;
t73 = t38 * t44;
t76 = (-t18 * t28 - t29 * t73) * mrSges(7,1) + (-t18 * t29 + t28 * t73) * mrSges(7,2);
t17 = t41 * t74 - t71 * t43;
t9 = t16 * t41 - t43 * t64;
t7 = t14 * t41 + t43 * t65;
t1 = [(-m(2) - m(3) - t85) * g(3) (-t85 * (pkin(2) * t73 + pkin(8) * t74) + (t94 * t42 + t97 * t44) * t38) * g(3) + (t93 * t13 + t80 * t14) * g(2) + (t93 * t15 + t80 * t16) * g(1) (t84 * t17 + t83 * t18) * g(3) + (t84 * t7 + t83 * t8) * g(2) + (t83 * t10 + t84 * t9) * g(1), t88 * (-g(1) * t9 - g(2) * t7 - g(3) * t17) (-(-t18 * t32 + t31 * t73) * mrSges(6,2) - t76 + t90 * (-t18 * t31 - t32 * t73)) * g(3) + (-(-t13 * t31 - t32 * t8) * mrSges(6,2) - t78 + t90 * (t13 * t32 - t31 * t8)) * g(2) + (-(-t10 * t32 - t15 * t31) * mrSges(6,2) - t77 + t90 * (-t10 * t31 + t15 * t32)) * g(1), -g(1) * t77 - g(2) * t78 - g(3) * t76];
taug  = t1(:);

% Calculate Gravitation load on the joints for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:13
% EndTime: 2019-03-08 20:45:15
% DurationCPUTime: 0.78s
% Computational Cost: add. (426->81), mult. (947->116), div. (0->0), fcn. (1088->12), ass. (0->40)
t87 = m(5) + m(6) + m(7);
t59 = m(4) + t87;
t81 = -m(7) * pkin(5) - mrSges(6,1);
t75 = mrSges(5,2) - m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) - mrSges(6,3) - mrSges(7,3);
t29 = qJ(5) + qJ(6);
t27 = sin(t29);
t28 = cos(t29);
t32 = sin(qJ(5));
t35 = cos(qJ(5));
t76 = -mrSges(5,1) - m(7) * (pkin(5) * t35 + pkin(4)) - t28 * mrSges(7,1) + t27 * mrSges(7,2) - m(6) * pkin(4) - t35 * mrSges(6,1) + t32 * mrSges(6,2);
t86 = -t27 * mrSges(7,1) - t35 * mrSges(6,2) - t28 * mrSges(7,2) + t81 * t32 - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t85 = t76 * t33 - t75 * t36 + mrSges(3,2) - mrSges(4,3);
t84 = pkin(2) * t59 + t87 * pkin(8) - t86;
t72 = -t59 * qJ(3) + t85;
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t30 = sin(pkin(11));
t62 = cos(pkin(6));
t56 = t30 * t62;
t61 = cos(pkin(11));
t17 = -t34 * t56 + t61 * t37;
t16 = t61 * t34 + t37 * t56;
t31 = sin(pkin(6));
t66 = t30 * t31;
t8 = t16 * t33 + t36 * t66;
t70 = (t17 * t28 - t27 * t8) * mrSges(7,1) + (-t17 * t27 - t28 * t8) * mrSges(7,2);
t49 = t62 * t61;
t14 = t30 * t34 - t37 * t49;
t55 = t31 * t61;
t10 = -t14 * t33 + t36 * t55;
t15 = t30 * t37 + t34 * t49;
t69 = (t10 * t27 + t15 * t28) * mrSges(7,1) + (t10 * t28 - t15 * t27) * mrSges(7,2);
t64 = t31 * t37;
t19 = -t33 * t64 + t62 * t36;
t65 = t31 * t34;
t68 = (-t19 * t27 + t28 * t65) * mrSges(7,1) + (-t19 * t28 - t27 * t65) * mrSges(7,2);
t63 = pkin(2) * t64 + qJ(3) * t65;
t1 = [(-m(2) - m(3) - t59) * g(3) (-m(4) * t63 - t87 * (pkin(8) * t64 + t63) + (t85 * t34 + t86 * t37) * t31) * g(3) + (t84 * t14 + t72 * t15) * g(2) + (t84 * t16 + t72 * t17) * g(1) (-g(1) * t16 - g(2) * t14 + g(3) * t64) * t59 (t75 * t19 + t76 * (-t62 * t33 - t36 * t64)) * g(3) + (t76 * (t14 * t36 + t33 * t55) - t75 * t10) * g(2) + (t75 * t8 + t76 * (t16 * t36 - t33 * t66)) * g(1) (-(-t19 * t35 - t32 * t65) * mrSges(6,2) - t68 + t81 * (-t19 * t32 + t35 * t65)) * g(3) + (-(t10 * t35 - t15 * t32) * mrSges(6,2) - t69 + t81 * (t10 * t32 + t15 * t35)) * g(2) + (-(-t17 * t32 - t35 * t8) * mrSges(6,2) - t70 + t81 * (t17 * t35 - t32 * t8)) * g(1), -g(1) * t70 - g(2) * t69 - g(3) * t68];
taug  = t1(:);

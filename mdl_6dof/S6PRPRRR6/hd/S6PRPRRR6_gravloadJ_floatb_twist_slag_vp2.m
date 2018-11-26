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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:06:36
% EndTime: 2018-11-23 15:06:36
% DurationCPUTime: 0.69s
% Computational Cost: add. (993->88), mult. (1136->117), div. (0->0), fcn. (1088->16), ass. (0->48)
t88 = m(5) + m(6);
t33 = qJ(5) + qJ(6);
t31 = sin(t33);
t32 = cos(t33);
t37 = sin(qJ(5));
t40 = cos(qJ(5));
t69 = m(4) + m(7) + t88;
t86 = pkin(2) * t69 + m(7) * (pkin(5) * t37 + pkin(8)) + t31 * mrSges(7,1) + t32 * mrSges(7,2) + t37 * mrSges(6,1) + t40 * mrSges(6,2) + mrSges(3,1) - mrSges(4,2) + mrSges(5,3) + t88 * pkin(8);
t84 = -mrSges(5,1) - m(7) * (pkin(5) * t40 + pkin(4)) - t32 * mrSges(7,1) + t31 * mrSges(7,2) - m(6) * pkin(4) - t40 * mrSges(6,1) + t37 * mrSges(6,2);
t87 = m(6) * pkin(9) - m(7) * (-pkin(10) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t85 = -m(7) * pkin(5) - mrSges(6,1);
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t79 = -t69 * qJ(3) + t84 * t38 + t87 * t41 + mrSges(3,2) - mrSges(4,3);
t34 = sin(pkin(11));
t70 = pkin(6) + qJ(2);
t61 = sin(t70);
t57 = t61 / 0.2e1;
t71 = pkin(6) - qJ(2);
t62 = sin(t71);
t52 = t57 - t62 / 0.2e1;
t42 = cos(qJ(2));
t72 = cos(pkin(11));
t64 = t72 * t42;
t18 = -t34 * t52 + t64;
t39 = sin(qJ(2));
t59 = cos(t71) / 0.2e1;
t63 = cos(t70);
t47 = t59 + t63 / 0.2e1;
t17 = t34 * t47 + t72 * t39;
t35 = sin(pkin(6));
t74 = t34 * t35;
t8 = t17 * t38 + t41 * t74;
t77 = (t18 * t32 - t31 * t8) * mrSges(7,1) + (-t18 * t31 - t32 * t8) * mrSges(7,2);
t14 = t34 * t39 - t72 * t47;
t65 = t35 * t72;
t10 = -t14 * t38 + t41 * t65;
t73 = t34 * t42;
t15 = t72 * t52 + t73;
t76 = (t10 * t31 + t15 * t32) * mrSges(7,1) + (t10 * t32 - t15 * t31) * mrSges(7,2);
t58 = t62 / 0.2e1;
t25 = t57 + t58;
t36 = cos(pkin(6));
t21 = -t25 * t38 + t36 * t41;
t26 = t59 - t63 / 0.2e1;
t75 = (-t21 * t31 + t26 * t32) * mrSges(7,1) + (-t21 * t32 - t26 * t31) * mrSges(7,2);
t53 = t58 - t61 / 0.2e1;
t1 = [(-m(2) - m(3) - t69) * g(3) (-t86 * t25 + t79 * t26) * g(3) + (t79 * (-t72 * t53 + t73) + t86 * t14) * g(2) + (t79 * (t34 * t53 + t64) + t86 * t17) * g(1) (-g(1) * t17 - g(2) * t14 + g(3) * t25) * t69 (-t87 * t21 + t84 * (-t25 * t41 - t36 * t38)) * g(3) + (t84 * (t14 * t41 + t38 * t65) + t87 * t10) * g(2) + (-t87 * t8 + t84 * (t17 * t41 - t38 * t74)) * g(1) (-(-t21 * t40 - t26 * t37) * mrSges(6,2) - t75 + t85 * (-t21 * t37 + t26 * t40)) * g(3) + (-(t10 * t40 - t15 * t37) * mrSges(6,2) - t76 + t85 * (t10 * t37 + t15 * t40)) * g(2) + (-(-t18 * t37 - t40 * t8) * mrSges(6,2) - t77 + t85 * (t18 * t40 - t37 * t8)) * g(1), -g(1) * t77 - g(2) * t76 - g(3) * t75];
taug  = t1(:);

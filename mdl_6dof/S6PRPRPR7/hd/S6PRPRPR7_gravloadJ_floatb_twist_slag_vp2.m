% Calculate Gravitation load on the joints for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2018-11-23 14:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:59:02
% EndTime: 2018-11-23 14:59:02
% DurationCPUTime: 0.69s
% Computational Cost: add. (804->88), mult. (953->112), div. (0->0), fcn. (890->14), ass. (0->54)
t76 = m(6) + m(7);
t62 = m(4) + m(5) + t76;
t72 = m(7) * pkin(9) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t77 = pkin(4) * t76 + t72;
t35 = sin(qJ(6));
t38 = cos(qJ(6));
t71 = t35 * mrSges(7,1) + t38 * mrSges(7,2) + qJ(5) * t76 - mrSges(5,2) + mrSges(6,3);
t78 = -m(7) * (pkin(5) + pkin(8)) - t38 * mrSges(7,1) + t35 * mrSges(7,2) - mrSges(3,1) - mrSges(6,1) - mrSges(5,3) + mrSges(4,2);
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t70 = -t62 * qJ(3) - t72 * t36 + t71 * t39 + mrSges(3,2) - mrSges(4,3);
t68 = pkin(4) * t36;
t32 = sin(pkin(10));
t33 = sin(pkin(6));
t67 = t32 * t33;
t40 = cos(qJ(2));
t66 = t32 * t40;
t65 = cos(pkin(10));
t64 = pkin(6) - qJ(2);
t63 = pkin(6) + qJ(2);
t37 = sin(qJ(2));
t53 = cos(t64) / 0.2e1;
t56 = cos(t63);
t43 = t53 + t56 / 0.2e1;
t14 = t32 * t37 - t65 * t43;
t11 = t14 * pkin(2);
t61 = -t14 * pkin(8) - t11;
t17 = t32 * t43 + t65 * t37;
t12 = t17 * pkin(2);
t60 = -t17 * pkin(8) - t12;
t54 = sin(t63);
t51 = t54 / 0.2e1;
t55 = sin(t64);
t52 = t55 / 0.2e1;
t27 = t51 + t52;
t26 = t27 * pkin(2);
t59 = t27 * pkin(8) + t26;
t58 = t33 * t65;
t57 = t65 * t40;
t49 = t52 - t54 / 0.2e1;
t48 = t51 - t55 / 0.2e1;
t34 = cos(pkin(6));
t28 = t53 - t56 / 0.2e1;
t22 = t28 * t68;
t20 = t27 * t39 + t34 * t36;
t19 = t32 * t49 + t57;
t18 = -t32 * t48 + t57;
t16 = -t65 * t49 + t66;
t15 = t65 * t48 + t66;
t8 = t19 * t68;
t7 = t16 * t68;
t5 = t14 * t39 + t36 * t58;
t3 = -t17 * t39 + t36 * t67;
t1 = [(-m(2) - m(3) - t62) * g(3) (-m(4) * t26 - m(5) * t59 - m(6) * (t22 + t59) - m(7) * (t22 + t26) + t70 * t28 + t78 * t27) * g(3) + (m(4) * t11 - m(5) * t61 - m(6) * (t61 + t7) - m(7) * (-t11 + t7) + t70 * t16 - t78 * t14) * g(2) + (m(4) * t12 - m(5) * t60 - m(6) * (t60 + t8) - m(7) * (-t12 + t8) + t70 * t19 - t78 * t17) * g(1) (-g(1) * t17 - g(2) * t14 + g(3) * t27) * t62 (-t71 * (-t27 * t36 + t34 * t39) + t77 * t20) * g(3) + (t71 * (-t14 * t36 + t39 * t58) - t77 * t5) * g(2) + (-t71 * (t17 * t36 + t39 * t67) + t77 * t3) * g(1), t76 * (-g(1) * t3 + g(2) * t5 - g(3) * t20) -g(1) * ((-t18 * t35 + t3 * t38) * mrSges(7,1) + (-t18 * t38 - t3 * t35) * mrSges(7,2)) - g(2) * ((-t15 * t35 - t38 * t5) * mrSges(7,1) + (-t15 * t38 + t35 * t5) * mrSges(7,2)) - g(3) * ((t20 * t38 - t28 * t35) * mrSges(7,1) + (-t20 * t35 - t28 * t38) * mrSges(7,2))];
taug  = t1(:);

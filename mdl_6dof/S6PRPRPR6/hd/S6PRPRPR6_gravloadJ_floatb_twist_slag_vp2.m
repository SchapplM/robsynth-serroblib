% Calculate Gravitation load on the joints for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 14:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:58:15
% EndTime: 2018-11-23 14:58:16
% DurationCPUTime: 0.59s
% Computational Cost: add. (860->76), mult. (987->99), div. (0->0), fcn. (929->16), ass. (0->47)
t75 = m(6) + m(7);
t63 = m(4) + m(5) + t75;
t29 = pkin(11) + qJ(6);
t27 = sin(t29);
t28 = cos(t29);
t30 = sin(pkin(11));
t33 = cos(pkin(11));
t76 = pkin(2) * t63 + m(7) * (pkin(5) * t30 + pkin(8)) + t27 * mrSges(7,1) + t28 * mrSges(7,2) + t30 * mrSges(6,1) + t33 * mrSges(6,2) + mrSges(3,1) - mrSges(4,2) + mrSges(5,3) + (m(5) + m(6)) * pkin(8);
t74 = mrSges(5,1) + m(7) * (pkin(5) * t33 + pkin(4)) + mrSges(7,1) * t28 - mrSges(7,2) * t27 + m(6) * pkin(4) + mrSges(6,1) * t33 - mrSges(6,2) * t30;
t77 = m(6) * qJ(5) - m(7) * (-pkin(9) - qJ(5)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t69 = -t63 * qJ(3) - t74 * t36 + t77 * t38 + mrSges(3,2) - mrSges(4,3);
t31 = sin(pkin(10));
t32 = sin(pkin(6));
t68 = t31 * t32;
t39 = cos(qJ(2));
t67 = t31 * t39;
t66 = cos(pkin(10));
t65 = pkin(6) - qJ(2);
t64 = pkin(6) + qJ(2);
t59 = t32 * t66;
t58 = t66 * t39;
t57 = cos(t64);
t56 = sin(t65);
t55 = sin(t64);
t54 = cos(t65) / 0.2e1;
t53 = t56 / 0.2e1;
t52 = t55 / 0.2e1;
t50 = t53 - t55 / 0.2e1;
t49 = t52 - t56 / 0.2e1;
t43 = t54 + t57 / 0.2e1;
t37 = sin(qJ(2));
t34 = cos(pkin(6));
t22 = t54 - t57 / 0.2e1;
t21 = t52 + t53;
t16 = -t21 * t36 + t34 * t38;
t15 = t21 * t38 + t34 * t36;
t13 = -t31 * t49 + t58;
t12 = t31 * t43 + t66 * t37;
t10 = t66 * t49 + t67;
t9 = t31 * t37 - t66 * t43;
t4 = -t9 * t36 + t38 * t59;
t3 = t36 * t59 + t9 * t38;
t2 = t12 * t36 + t38 * t68;
t1 = -t12 * t38 + t36 * t68;
t5 = [(-m(2) - m(3) - t63) * g(3) (-t76 * t21 + t69 * t22) * g(3) + (t69 * (-t66 * t50 + t67) + t76 * t9) * g(2) + (t69 * (t31 * t50 + t58) + t76 * t12) * g(1) (-g(1) * t12 - g(2) * t9 + g(3) * t21) * t63 (t74 * t15 - t16 * t77) * g(3) + (-t74 * t3 + t4 * t77) * g(2) + (t74 * t1 - t2 * t77) * g(1), t75 * (-g(1) * t1 + g(2) * t3 - g(3) * t15) -g(1) * ((t13 * t28 - t2 * t27) * mrSges(7,1) + (-t13 * t27 - t2 * t28) * mrSges(7,2)) - g(2) * ((t10 * t28 + t27 * t4) * mrSges(7,1) + (-t10 * t27 + t28 * t4) * mrSges(7,2)) - g(3) * ((-t16 * t27 + t22 * t28) * mrSges(7,1) + (-t16 * t28 - t22 * t27) * mrSges(7,2))];
taug  = t5(:);

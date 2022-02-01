% Calculate Gravitation load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:51
% EndTime: 2022-01-23 09:31:53
% DurationCPUTime: 0.50s
% Computational Cost: add. (232->83), mult. (329->99), div. (0->0), fcn. (308->8), ass. (0->48)
t81 = mrSges(5,1) + mrSges(6,1);
t80 = mrSges(5,2) + mrSges(6,2);
t37 = cos(qJ(3));
t29 = t37 * pkin(3);
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t58 = pkin(2) * t34 + pkin(1);
t65 = pkin(7) + pkin(6);
t79 = m(4) * t58 + m(5) * (t29 * t34 + t58) + mrSges(2,1) + t34 * mrSges(3,1) + (m(4) * pkin(6) + m(5) * t65 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3)) * t33;
t77 = m(5) * pkin(3) + mrSges(4,1);
t32 = qJ(3) + qJ(4);
t24 = sin(t32);
t25 = cos(t32);
t73 = mrSges(5,1) * t24 + t80 * t25;
t38 = cos(qJ(1));
t48 = t38 * t24;
t36 = sin(qJ(1));
t53 = t36 * t25;
t12 = -t34 * t48 + t53;
t47 = t38 * t25;
t54 = t36 * t24;
t13 = t34 * t47 + t54;
t72 = -t81 * t12 + t80 * t13;
t10 = t34 * t54 + t47;
t11 = -t34 * t53 + t48;
t71 = t81 * t10 - t80 * t11;
t70 = -m(3) - m(4) - m(6);
t35 = sin(qJ(3));
t59 = t35 * pkin(3);
t69 = -m(5) * (qJ(2) + t59) + mrSges(2,2) - mrSges(3,3);
t66 = m(6) * pkin(4);
t60 = g(3) * t33;
t19 = pkin(4) * t24 + t59;
t55 = t36 * t19;
t51 = t36 * t35;
t50 = t36 * t37;
t49 = t38 * t19;
t45 = t38 * t35;
t44 = t38 * t37;
t20 = pkin(4) * t25 + t29;
t26 = t36 * qJ(2);
t43 = t38 * pkin(1) + t26;
t39 = (pkin(2) + t20) * t34 - (-qJ(5) - t65) * t33;
t16 = -t34 * t45 + t50;
t14 = t34 * t51 + t44;
t17 = t34 * t44 + t51;
t15 = -t34 * t50 + t45;
t1 = [(-m(3) * t43 - m(4) * t26 - t17 * mrSges(4,1) - t16 * mrSges(4,2) - m(6) * (t43 + t55) - t81 * t13 - t80 * t12 + t69 * t36) * g(2) + (-m(6) * t49 - t15 * mrSges(4,1) - t14 * mrSges(4,2) - t80 * t10 - t81 * t11 + (m(3) * pkin(1) - m(6) * (-pkin(1) - t39) + t79) * t36) * g(1) + ((-m(6) * t39 - t79) * g(2) + (t70 * qJ(2) + t69) * g(1)) * t38, (-g(1) * t36 + g(2) * t38) * (m(5) - t70), (m(5) * t59 + m(6) * t19 + mrSges(4,1) * t35 + mrSges(6,1) * t24 + mrSges(4,2) * t37 + t73) * t60 + (-t15 * mrSges(4,2) - m(6) * (-t38 * t20 - t34 * t55) + t77 * t14 + t71) * g(2) + (t17 * mrSges(4,2) - m(6) * (t36 * t20 - t34 * t49) - t77 * t16 + t72) * g(1), (-(-mrSges(6,1) - t66) * t24 + t73) * t60 + (t10 * t66 + t71) * g(2) + (-t12 * t66 + t72) * g(1), (g(3) * t34 + (-g(1) * t38 - g(2) * t36) * t33) * m(6)];
taug = t1(:);

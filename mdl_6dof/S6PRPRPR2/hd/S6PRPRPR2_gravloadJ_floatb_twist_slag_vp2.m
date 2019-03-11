% Calculate Gravitation load on the joints for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:35
% EndTime: 2019-03-08 19:29:37
% DurationCPUTime: 0.87s
% Computational Cost: add. (572->87), mult. (1372->130), div. (0->0), fcn. (1681->14), ass. (0->51)
t34 = pkin(12) + qJ(6);
t32 = sin(t34);
t33 = cos(t34);
t35 = sin(pkin(12));
t38 = cos(pkin(12));
t90 = m(7) * (pkin(5) * t38 + pkin(4)) + t33 * mrSges(7,1) - t32 * mrSges(7,2) + m(6) * pkin(4) + t38 * mrSges(6,1) - t35 * mrSges(6,2) + mrSges(5,1);
t89 = m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t91 = m(6) + m(7);
t96 = -m(5) - t91;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t94 = t89 * t42 - t90 * t44 - mrSges(4,1);
t92 = m(7) * (pkin(5) * t35 + pkin(8)) + t32 * mrSges(7,1) + t33 * mrSges(7,2) + t35 * mrSges(6,1) + t38 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3) + (m(5) + m(6)) * pkin(8);
t37 = sin(pkin(10));
t39 = cos(pkin(10));
t43 = sin(qJ(2));
t40 = cos(pkin(6));
t85 = cos(qJ(2));
t72 = t40 * t85;
t97 = -t37 * t43 + t39 * t72;
t36 = sin(pkin(11));
t76 = cos(pkin(11));
t21 = -t85 * t36 - t43 * t76;
t88 = m(4) - t96;
t52 = -t43 * t36 + t76 * t85;
t77 = t21 * t40;
t12 = t37 * t77 + t39 * t52;
t7 = -t37 * t52 + t39 * t77;
t81 = t40 * t43;
t75 = sin(pkin(6));
t68 = t42 * t75;
t66 = t43 * t75;
t65 = t44 * t75;
t63 = t97 * pkin(2);
t62 = t85 * t75;
t58 = t76 * t75;
t28 = pkin(2) * t62;
t54 = -t37 * t72 - t39 * t43;
t51 = t54 * pkin(2);
t50 = t40 * t52;
t19 = t36 * t62 + t43 * t58;
t18 = t36 * t66 - t58 * t85;
t14 = t19 * t44 + t40 * t42;
t13 = t19 * t42 - t40 * t44;
t11 = t21 * t39 - t37 * t50;
t8 = t37 * t21 + t39 * t50;
t4 = t12 * t44 + t37 * t68;
t3 = t12 * t42 - t37 * t65;
t2 = -t39 * t68 - t44 * t7;
t1 = t39 * t65 - t42 * t7;
t5 = [(-m(2) - m(3) - t88) * g(3) (-m(4) * t28 - mrSges(3,1) * t62 + mrSges(3,2) * t66 + t96 * (-t18 * pkin(3) + t28) - t92 * t19 - t94 * t18) * g(3) + (-t97 * mrSges(3,1) - (-t37 * t85 - t39 * t81) * mrSges(3,2) - m(4) * t63 + t96 * (t8 * pkin(3) + t63) + t94 * t8 + t92 * t7) * g(2) + (-t54 * mrSges(3,1) - (t37 * t81 - t39 * t85) * mrSges(3,2) - m(4) * t51 + t96 * (t11 * pkin(3) + t51) + t94 * t11 - t92 * t12) * g(1) (-g(3) * t40 + (-g(1) * t37 + g(2) * t39) * t75) * t88 (t90 * t13 + t89 * t14) * g(3) + (t90 * t1 + t89 * t2) * g(2) + (t90 * t3 + t89 * t4) * g(1), t91 * (-g(1) * t3 - g(2) * t1 - g(3) * t13) -g(1) * ((-t11 * t33 - t32 * t4) * mrSges(7,1) + (t11 * t32 - t33 * t4) * mrSges(7,2)) - g(2) * ((-t2 * t32 - t33 * t8) * mrSges(7,1) + (-t2 * t33 + t32 * t8) * mrSges(7,2)) - g(3) * ((-t14 * t32 + t18 * t33) * mrSges(7,1) + (-t14 * t33 - t18 * t32) * mrSges(7,2))];
taug  = t5(:);

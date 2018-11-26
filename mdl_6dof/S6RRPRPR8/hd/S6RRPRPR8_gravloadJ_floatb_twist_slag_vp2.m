% Calculate Gravitation load on the joints for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:05:15
% EndTime: 2018-11-23 17:05:16
% DurationCPUTime: 1.05s
% Computational Cost: add. (529->125), mult. (745->141), div. (0->0), fcn. (750->10), ass. (0->66)
t96 = mrSges(5,1) + mrSges(6,1);
t95 = mrSges(5,2) - mrSges(6,3);
t109 = -m(4) * qJ(3) - mrSges(4,3) + mrSges(7,3);
t33 = -pkin(8) - qJ(3);
t108 = -m(7) * (-pkin(9) - t33) + t109;
t107 = mrSges(5,3) + mrSges(6,2);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t102 = g(1) * t39 + g(2) * t36;
t35 = sin(qJ(2));
t106 = t102 * t35;
t30 = pkin(10) + qJ(4);
t25 = sin(t30);
t26 = cos(t30);
t34 = sin(qJ(6));
t37 = cos(qJ(6));
t51 = t25 * t34 + t26 * t37;
t52 = t25 * t37 - t26 * t34;
t105 = t51 * mrSges(7,1) + t52 * mrSges(7,2) - t95 * t25 + t96 * t26;
t38 = cos(qJ(2));
t77 = t36 * t25;
t10 = t26 * t39 + t38 * t77;
t74 = t39 * t25;
t76 = t38 * t26;
t11 = t36 * t76 - t74;
t53 = t10 * t34 + t11 * t37;
t99 = t10 * t37 - t11 * t34;
t104 = t99 * mrSges(7,1) - t53 * mrSges(7,2);
t31 = sin(pkin(10));
t32 = cos(pkin(10));
t47 = m(4) * pkin(2) + t32 * mrSges(4,1) - t31 * mrSges(4,2);
t98 = t47 * t38;
t103 = t108 * t35 - t98;
t101 = m(3) + m(4);
t100 = m(6) + m(7);
t72 = qJ(5) * t25;
t94 = pkin(4) * t76 + t38 * t72;
t93 = -m(5) - t100;
t92 = -mrSges(4,1) * t31 - mrSges(4,2) * t32 + mrSges(2,2) - mrSges(3,3);
t58 = t38 * mrSges(3,1) - t35 * mrSges(3,2);
t91 = t107 * t35 + t58;
t90 = m(7) * pkin(5) + t96;
t87 = pkin(3) * t31;
t86 = pkin(5) * t26;
t82 = t25 * t35;
t81 = t33 * t35;
t78 = t35 * t39;
t24 = pkin(3) * t32 + pkin(2);
t17 = t38 * t24;
t75 = t38 * t39;
t73 = t39 * pkin(1) + t36 * pkin(7);
t71 = t33 * t78;
t28 = t39 * pkin(7);
t70 = t36 * t81 + t39 * t87 + t28;
t63 = t17 - t81;
t62 = t24 * t75 + t36 * t87 + t73;
t12 = -t36 * t26 + t38 * t74;
t13 = t26 * t75 + t77;
t1 = t12 * t37 - t13 * t34;
t2 = t12 * t34 + t13 * t37;
t60 = mrSges(7,1) * t1 - mrSges(7,2) * t2;
t59 = (mrSges(7,1) * t52 - mrSges(7,2) * t51) * t35;
t50 = -pkin(4) * t26 - t24 - t72;
t46 = t13 * pkin(4) + t12 * qJ(5) + t62;
t15 = t35 * t26 * qJ(5);
t3 = [(-m(5) * (t62 - t71) - m(6) * (t46 - t71) - m(7) * t46 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t107 * t78 - t101 * t73 - t90 * t13 + t95 * t12 + t92 * t36 + (-mrSges(2,1) - t58 + t103) * t39) * g(2) + (-m(5) * t70 + t53 * mrSges(7,1) + t99 * mrSges(7,2) - t100 * (-t11 * pkin(4) - qJ(5) * t10 + t70) - t101 * t28 + t90 * t11 - t95 * t10 + t92 * t39 + (t98 + mrSges(2,1) + t101 * pkin(1) + t93 * (-pkin(1) - t17) + (-m(7) * pkin(9) - t109) * t35 + t91) * t36) * g(1) (-m(5) * t63 - m(6) * (t63 + t94) - m(7) * (t17 + t94) - t91 + t103) * g(3) + ((-m(7) * t86 - t105) * g(3) + t102 * (mrSges(3,2) - t107 + (m(5) + m(6)) * t33 + t108)) * t38 + (mrSges(3,1) + t47 + m(5) * t24 - m(6) * t50 - m(7) * (t50 - t86) + t105) * t106 (t38 * g(3) - t106) * (m(4) - t93) (-m(6) * t15 - m(7) * (t15 + (-pkin(4) - pkin(5)) * t82) + t59 + (t95 * t26 + (m(6) * pkin(4) + t96) * t25) * t35) * g(3) + (-t100 * (-t10 * pkin(4) + qJ(5) * t11) + t95 * t11 + t90 * t10 + t104) * g(2) + (t60 - t100 * (-t12 * pkin(4) + qJ(5) * t13) + t95 * t13 + t90 * t12) * g(1), t100 * (-g(1) * t12 - g(2) * t10 - g(3) * t82) -g(1) * t60 - g(2) * t104 - g(3) * t59];
taug  = t3(:);

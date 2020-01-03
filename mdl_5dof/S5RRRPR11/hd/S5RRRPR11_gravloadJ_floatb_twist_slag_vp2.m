% Calculate Gravitation load on the joints for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:26
% EndTime: 2019-12-31 21:32:29
% DurationCPUTime: 0.92s
% Computational Cost: add. (263->107), mult. (613->130), div. (0->0), fcn. (636->8), ass. (0->52)
t80 = mrSges(4,1) + mrSges(5,1);
t78 = mrSges(4,2) - mrSges(5,3);
t87 = mrSges(5,2) + mrSges(4,3);
t86 = m(6) * pkin(8) + mrSges(6,3) - t87;
t29 = sin(qJ(3));
t33 = cos(qJ(3));
t28 = sin(qJ(5));
t32 = cos(qJ(5));
t43 = t28 * t29 + t32 * t33;
t44 = t28 * t33 - t29 * t32;
t85 = t43 * mrSges(6,1) - t44 * mrSges(6,2) - t78 * t29 + t80 * t33;
t35 = cos(qJ(1));
t31 = sin(qJ(1));
t34 = cos(qJ(2));
t65 = t31 * t34;
t10 = t29 * t65 + t33 * t35;
t62 = t35 * t29;
t64 = t33 * t34;
t11 = t31 * t64 - t62;
t45 = t10 * t28 + t11 * t32;
t82 = t10 * t32 - t11 * t28;
t84 = t82 * mrSges(6,1) - t45 * mrSges(6,2);
t83 = m(5) + m(6);
t79 = mrSges(2,2) - mrSges(3,3);
t30 = sin(qJ(2));
t49 = t34 * mrSges(3,1) - t30 * mrSges(3,2);
t77 = t30 * t87 + t49;
t75 = pkin(7) * (-m(4) - t83);
t59 = qJ(4) * t29;
t41 = -pkin(3) * t33 - pkin(2) - t59;
t71 = pkin(4) * t33;
t74 = (mrSges(3,2) + t86) * t34 + (-m(6) * (t41 - t71) - m(5) * t41 + m(4) * pkin(2) + mrSges(3,1) + t85) * t30;
t73 = m(6) * pkin(4) + t80;
t23 = t30 * pkin(7);
t25 = t34 * pkin(2);
t69 = t29 * t30;
t66 = t30 * t35;
t63 = t34 * t35;
t61 = t25 + t23;
t60 = pkin(1) * t35 + pkin(6) * t31;
t58 = -pkin(1) - t25;
t54 = pkin(3) * t64 + t34 * t59 + t61;
t53 = pkin(2) * t63 + pkin(7) * t66 + t60;
t12 = -t31 * t33 + t34 * t62;
t13 = t29 * t31 + t33 * t63;
t1 = t12 * t32 - t13 * t28;
t2 = t12 * t28 + t13 * t32;
t51 = mrSges(6,1) * t1 - mrSges(6,2) * t2;
t50 = (-mrSges(6,1) * t44 - mrSges(6,2) * t43) * t30;
t26 = t35 * pkin(6);
t14 = t33 * t30 * qJ(4);
t3 = [(-m(3) * t60 - m(4) * t53 - t2 * mrSges(6,1) - t1 * mrSges(6,2) - t83 * (pkin(3) * t13 + qJ(4) * t12 + t53) + (-mrSges(2,1) - t49) * t35 + t79 * t31 - t73 * t13 + t78 * t12 + t86 * t66) * g(2) + (t45 * mrSges(6,1) + t82 * mrSges(6,2) - t83 * (-pkin(3) * t11 - qJ(4) * t10 + t26) + t79 * t35 + (-m(3) - m(4)) * t26 + t73 * t11 - t78 * t10 + (mrSges(2,1) + m(3) * pkin(1) - m(6) * t58 - (m(6) * (-pkin(7) + pkin(8)) + mrSges(6,3)) * t30 + (-m(4) - m(5)) * (t58 - t23) + t77) * t31) * g(1), (-m(4) * t61 - m(5) * t54 - m(6) * (-pkin(8) * t30 + t54) + t30 * mrSges(6,3) + (-m(6) * t71 - t85) * t34 - t77) * g(3) + (t31 * t74 + t65 * t75) * g(2) + (t35 * t74 + t63 * t75) * g(1), (-m(5) * t14 - m(6) * (t14 + (-pkin(3) - pkin(4)) * t69) + t50 + (t78 * t33 + (m(5) * pkin(3) + t80) * t29) * t30) * g(3) + (-t83 * (-pkin(3) * t10 + qJ(4) * t11) + t78 * t11 + t73 * t10 + t84) * g(2) + (t51 - t83 * (-pkin(3) * t12 + qJ(4) * t13) + t78 * t13 + t73 * t12) * g(1), t83 * (-g(1) * t12 - g(2) * t10 - g(3) * t69), -g(1) * t51 - g(2) * t84 - g(3) * t50];
taug = t3(:);

% Calculate Gravitation load on the joints for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:40
% EndTime: 2019-12-31 21:50:41
% DurationCPUTime: 0.39s
% Computational Cost: add. (323->76), mult. (272->83), div. (0->0), fcn. (213->8), ass. (0->41)
t85 = -mrSges(5,1) - mrSges(6,1);
t36 = qJ(3) + qJ(4);
t30 = sin(t36);
t87 = (mrSges(5,2) - mrSges(6,3)) * t30;
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t82 = t40 * mrSges(4,1) - t38 * mrSges(4,2);
t86 = -mrSges(3,1) - t82;
t32 = cos(t36);
t80 = t85 * t32 + t87;
t84 = -mrSges(6,2) - mrSges(4,3) - mrSges(5,3) + mrSges(3,2);
t81 = m(6) * qJ(5) * t32;
t14 = t30 * qJ(5);
t64 = t32 * pkin(4) + t14;
t37 = qJ(1) + qJ(2);
t31 = sin(t37);
t33 = cos(t37);
t79 = g(1) * t33 + g(2) * t31;
t34 = t40 * pkin(3);
t29 = t34 + pkin(2);
t78 = t84 * t33 + (-m(6) * (-t29 - t64) - t80 - t86) * t31;
t69 = t32 * t33;
t77 = t84 * t31 + t85 * t69 + (t87 + t86) * t33;
t76 = pkin(3) * t38;
t39 = sin(qJ(1));
t73 = t39 * pkin(1);
t41 = cos(qJ(1));
t35 = t41 * pkin(1);
t72 = mrSges(5,2) * t32;
t42 = -pkin(8) - pkin(7);
t68 = t33 * t42;
t63 = t33 * pkin(2) + t31 * pkin(7);
t60 = (-mrSges(6,3) * t32 - t81) * t31;
t59 = -mrSges(6,3) * t69 - t33 * t81;
t58 = -t31 * pkin(2) + t33 * pkin(7);
t57 = t33 * t29 - t31 * t42;
t52 = -t31 * t29 - t68;
t51 = pkin(4) * t69 + t33 * t14 + t57;
t45 = t72 + (m(6) * pkin(4) - t85) * t30;
t44 = m(6) * (-pkin(4) * t30 - t76) - t30 * mrSges(6,1);
t1 = [(-t41 * mrSges(2,1) + t39 * mrSges(2,2) - m(3) * t35 - m(4) * (t35 + t63) - m(5) * (t35 + t57) - m(6) * (t35 + t51) + t77) * g(2) + (t39 * mrSges(2,1) + t41 * mrSges(2,2) + m(3) * t73 - m(4) * (t58 - t73) - m(5) * (t52 - t73) - m(6) * (-t68 - t73) + t78) * g(1), (-m(4) * t63 - m(5) * t57 - m(6) * t51 + t77) * g(2) + (-m(4) * t58 - m(5) * t52 + m(6) * t68 + t78) * g(1), -g(1) * (t44 * t33 - t59) - g(2) * (t44 * t31 - t60) + (-m(5) * t34 - m(6) * (t34 + t64) + t80 - t82) * g(3) + (m(5) * t76 + mrSges(4,1) * t38 + mrSges(5,1) * t30 + mrSges(4,2) * t40 + t72) * t79, (-m(6) * t64 + t80) * g(3) + (t45 * t31 + t60) * g(2) + (t45 * t33 + t59) * g(1), (g(3) * t32 - t79 * t30) * m(6)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:53
% EndTime: 2019-12-31 20:16:55
% DurationCPUTime: 0.44s
% Computational Cost: add. (296->76), mult. (294->84), div. (0->0), fcn. (245->10), ass. (0->42)
t77 = m(5) + m(6);
t28 = qJ(2) + pkin(9);
t25 = qJ(4) + t28;
t20 = sin(t25);
t21 = cos(t25);
t33 = cos(qJ(5));
t60 = t33 * mrSges(6,1);
t36 = mrSges(5,2) * t21 + (m(6) * pkin(4) + mrSges(5,1) + t60) * t20;
t30 = sin(qJ(5));
t63 = t30 * mrSges(6,2);
t84 = -t20 * t63 + t21 * (-m(6) * pkin(8) - mrSges(6,3));
t75 = t21 * pkin(4) + t20 * pkin(8);
t83 = m(6) * t75;
t82 = -t21 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t20;
t23 = sin(t28);
t24 = cos(t28);
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t68 = t31 * pkin(2);
t81 = m(4) * t68 + mrSges(3,1) * t31 + mrSges(4,1) * t23 + mrSges(3,2) * t34 + mrSges(4,2) * t24 - t77 * (-pkin(3) * t23 - t68) + t36;
t80 = -t34 * mrSges(3,1) - t24 * mrSges(4,1) + t31 * mrSges(3,2) + t23 * mrSges(4,2);
t74 = -(t60 - t63) * t21 + t82;
t29 = -qJ(3) - pkin(6);
t72 = -m(3) * pkin(6) + m(4) * t29 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t26 = t34 * pkin(2);
t71 = mrSges(2,1) + m(4) * (t26 + pkin(1)) + m(3) * pkin(1) - t80 - t82;
t32 = sin(qJ(1));
t62 = t32 * t30;
t61 = t32 * t33;
t35 = cos(qJ(1));
t59 = t35 * t30;
t58 = t35 * t33;
t56 = pkin(3) * t24 + t26;
t49 = t84 * t32;
t48 = t84 * t35;
t27 = -pkin(7) + t29;
t8 = pkin(1) + t56;
t4 = t21 * t58 + t62;
t3 = -t21 * t59 + t61;
t2 = -t21 * t61 + t59;
t1 = t21 * t62 + t58;
t5 = [(-t4 * mrSges(6,1) - t3 * mrSges(6,2) - t77 * (-t32 * t27 + t35 * t8) + t72 * t32 + (-t71 - t83) * t35) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (t77 * t27 + t72) * t35 + (m(5) * t8 - m(6) * (-t75 - t8) + t71) * t32) * g(1), (-m(4) * t26 - m(5) * t56 - m(6) * (t56 + t75) + t74 + t80) * g(3) + (t81 * t32 + t49) * g(2) + (t81 * t35 + t48) * g(1), (-g(1) * t32 + g(2) * t35) * (m(4) + t77), (t74 - t83) * g(3) + (t36 * t32 + t49) * g(2) + (t36 * t35 + t48) * g(1), -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t30 - mrSges(6,2) * t33) * t20];
taug = t5(:);

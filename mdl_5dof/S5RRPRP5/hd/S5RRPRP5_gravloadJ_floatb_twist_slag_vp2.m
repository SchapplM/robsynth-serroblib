% Calculate Gravitation load on the joints for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:01
% EndTime: 2019-12-31 19:54:02
% DurationCPUTime: 0.38s
% Computational Cost: add. (257->66), mult. (247->68), div. (0->0), fcn. (190->8), ass. (0->33)
t69 = -mrSges(5,1) - mrSges(6,1);
t22 = qJ(2) + pkin(8);
t17 = sin(t22);
t18 = cos(t22);
t24 = sin(qJ(2));
t26 = cos(qJ(2));
t19 = qJ(4) + t22;
t14 = sin(t19);
t15 = cos(t19);
t59 = t69 * t15 + (mrSges(5,2) - mrSges(6,3)) * t14;
t68 = t26 * mrSges(3,1) + t18 * mrSges(4,1) - t24 * mrSges(3,2) - t17 * mrSges(4,2) - t59;
t61 = t15 * pkin(4) + t14 * qJ(5);
t67 = m(6) * t61;
t65 = (-m(6) * qJ(5) - mrSges(6,3)) * t15;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t64 = g(1) * t27 + g(2) * t25;
t63 = m(5) + m(6);
t20 = t26 * pkin(2);
t57 = mrSges(2,1) + m(4) * (t20 + pkin(1)) + m(3) * pkin(1) + t68;
t23 = -qJ(3) - pkin(6);
t56 = -m(3) * pkin(6) + m(4) * t23 + mrSges(2,2) - mrSges(6,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t53 = t24 * pkin(2);
t51 = mrSges(5,2) * t15;
t48 = pkin(3) * t18 + t20;
t45 = t65 * t25;
t44 = t65 * t27;
t3 = -pkin(3) * t17 - t53;
t30 = m(6) * (-pkin(4) * t14 + t3) - t14 * mrSges(6,1);
t28 = t51 + (m(6) * pkin(4) - t69) * t14;
t21 = -pkin(7) + t23;
t2 = pkin(1) + t48;
t1 = [(-t63 * (t27 * t2 - t25 * t21) + (-t57 - t67) * t27 + t56 * t25) * g(2) + ((t63 * t21 + t56) * t27 + (m(5) * t2 - m(6) * (-t2 - t61) + t57) * t25) * g(1), -g(1) * (t30 * t27 - t44) - g(2) * (t30 * t25 - t45) + (-m(4) * t20 - m(5) * t48 - m(6) * (t48 + t61) - t68) * g(3) + t64 * (m(4) * t53 - m(5) * t3 + mrSges(3,1) * t24 + mrSges(4,1) * t17 + mrSges(5,1) * t14 + mrSges(3,2) * t26 + mrSges(4,2) * t18 + t51), (-g(1) * t25 + g(2) * t27) * (m(4) + t63), (t59 - t67) * g(3) + (t28 * t25 + t45) * g(2) + (t28 * t27 + t44) * g(1), (g(3) * t15 - t14 * t64) * m(6)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:21
% EndTime: 2019-12-31 20:03:23
% DurationCPUTime: 0.61s
% Computational Cost: add. (180->77), mult. (395->88), div. (0->0), fcn. (368->6), ass. (0->39)
t58 = m(6) * pkin(4);
t72 = mrSges(5,1) + mrSges(6,1);
t71 = mrSges(5,2) + mrSges(6,2);
t24 = sin(qJ(2));
t27 = cos(qJ(2));
t70 = (-mrSges(3,1) - mrSges(4,1)) * t27 + (mrSges(3,2) - mrSges(4,3)) * t24;
t23 = sin(qJ(4));
t51 = t27 * t23;
t26 = cos(qJ(4));
t16 = t26 * pkin(4) + pkin(3);
t54 = -pkin(2) - t16;
t57 = -pkin(2) - pkin(3);
t61 = m(4) + m(5) + m(6);
t69 = -t51 * t58 + (m(4) * pkin(2) - m(5) * t57 - m(6) * t54 + mrSges(4,1)) * t24 + (-qJ(3) * t61 - mrSges(4,3)) * t27;
t68 = mrSges(2,1) - t70;
t67 = pkin(7) * m(5) - (-qJ(5) - pkin(7)) * m(6) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t28 = cos(qJ(1));
t50 = t27 * t28;
t52 = t24 * t26;
t3 = t23 * t50 - t28 * t52;
t53 = t24 * t23;
t33 = t27 * t26 + t53;
t4 = t33 * t28;
t66 = -t72 * t3 - t71 * t4;
t25 = sin(qJ(1));
t63 = t51 - t52;
t1 = t63 * t25;
t2 = t33 * t25;
t65 = -t72 * t1 - t71 * t2;
t64 = -t72 * t33 + t71 * t63;
t62 = g(1) * t28 + g(2) * t25;
t19 = t27 * pkin(2);
t17 = t24 * qJ(3);
t49 = t19 + t17;
t48 = t28 * pkin(1) + t25 * pkin(6);
t44 = -pkin(1) - t17;
t32 = pkin(4) * t53 + t27 * t16;
t31 = t63 * t58;
t5 = [(-m(5) * pkin(3) * t50 - m(3) * t48 - t61 * (pkin(2) * t50 + t28 * t17 + t48) - t72 * t4 + t71 * t3 + (-m(6) * t32 - t68) * t28 + t67 * t25) * g(2) + (t72 * t2 - t71 * t1 + (m(3) * pkin(1) - m(4) * (t44 - t19) - m(5) * (t57 * t27 + t44) - m(6) * (-pkin(1) + t54 * t27 + (-pkin(4) * t23 - qJ(3)) * t24) + t68) * t25 + ((-m(3) - t61) * pkin(6) + t67) * t28) * g(1), t62 * (mrSges(3,1) * t24 + mrSges(3,2) * t27) + (t69 * t25 + t65) * g(2) + (t69 * t28 + t66) * g(1) + (-m(4) * t49 - m(5) * (t27 * pkin(3) + t49) - m(6) * (t32 + t49) + t64 + t70) * g(3), (t27 * g(3) - t62 * t24) * t61, (t33 * t58 - t64) * g(3) + (t25 * t31 - t65) * g(2) + (t28 * t31 - t66) * g(1), (g(1) * t25 - g(2) * t28) * m(6)];
taug = t5(:);

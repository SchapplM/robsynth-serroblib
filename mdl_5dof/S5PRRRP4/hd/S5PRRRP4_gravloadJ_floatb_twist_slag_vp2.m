% Calculate Gravitation load on the joints for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:37
% EndTime: 2019-12-05 16:45:38
% DurationCPUTime: 0.46s
% Computational Cost: add. (239->61), mult. (324->74), div. (0->0), fcn. (285->8), ass. (0->39)
t86 = mrSges(5,1) + mrSges(6,1);
t28 = qJ(2) + qJ(3);
t25 = sin(t28);
t26 = cos(t28);
t31 = sin(qJ(4));
t54 = t31 * mrSges(5,2);
t81 = -mrSges(6,2) - mrSges(5,3);
t82 = -m(5) - m(6);
t85 = -t25 * t54 + t26 * (pkin(7) * t82 + t81);
t33 = cos(qJ(4));
t84 = t31 * mrSges(6,3) + t86 * t33;
t78 = pkin(4) * t33 + qJ(5) * t31;
t79 = (-m(6) * (-pkin(3) - t78) + t84) * t25;
t65 = pkin(3) * t25;
t32 = sin(qJ(2));
t66 = pkin(2) * t32;
t76 = m(6) * t66 - m(5) * (-t65 - t66) + t79;
t29 = sin(pkin(8));
t30 = cos(pkin(8));
t75 = g(1) * t30 + g(2) * t29;
t74 = (t54 - mrSges(4,1) - t84) * t26 + (mrSges(4,2) + t81) * t25;
t73 = m(6) * pkin(4) + t86;
t72 = t85 * t29;
t71 = t85 * t30;
t70 = m(5) * t65 + t79;
t69 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t61 = g(3) * t25;
t34 = cos(qJ(2));
t27 = t34 * pkin(2);
t58 = t29 * t31;
t57 = t29 * t33;
t56 = t30 * t31;
t55 = t30 * t33;
t50 = t26 * pkin(3) + t25 * pkin(7);
t42 = t78 * t26 + t50;
t40 = mrSges(4,1) * t25 + mrSges(4,2) * t26;
t3 = t26 * t56 - t57;
t1 = t26 * t58 + t55;
t2 = [(-m(2) - m(3) - m(4) + t82) * g(3), (t76 * t29 + t72) * g(2) + (t76 * t30 + t71) * g(1) + (-t34 * mrSges(3,1) + t32 * mrSges(3,2) - m(4) * t27 - m(5) * (t27 + t50) - m(6) * (t27 + t42) + t74) * g(3) + (m(4) * t66 + mrSges(3,1) * t32 + mrSges(3,2) * t34 + t40) * t75, t75 * t40 + (t70 * t29 + t72) * g(2) + (t70 * t30 + t71) * g(1) + (-m(5) * t50 - m(6) * t42 + t74) * g(3), (t73 * t31 + t69 * t33) * t61 + (t69 * (t26 * t57 - t56) + t73 * t1) * g(2) + (t69 * (t26 * t55 + t58) + t73 * t3) * g(1), (-g(1) * t3 - g(2) * t1 - t31 * t61) * m(6)];
taug = t2(:);

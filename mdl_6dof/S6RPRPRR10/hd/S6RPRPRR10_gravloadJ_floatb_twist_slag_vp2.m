% Calculate Gravitation load on the joints for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:39
% EndTime: 2019-03-09 04:07:41
% DurationCPUTime: 0.64s
% Computational Cost: add. (323->83), mult. (429->94), div. (0->0), fcn. (379->10), ass. (0->46)
t30 = -pkin(8) - qJ(4);
t78 = mrSges(4,2) + m(6) * t30 + m(7) * (-pkin(9) + t30) - mrSges(6,3) - mrSges(7,3) - m(5) * qJ(4) - mrSges(5,3);
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t77 = t31 * mrSges(4,1) + t78 * t33;
t29 = cos(pkin(10));
t17 = t29 * pkin(4) + pkin(3);
t27 = pkin(10) + qJ(5);
t19 = cos(t27);
t28 = sin(pkin(10));
t75 = -m(6) * t17 - m(7) * (pkin(5) * t19 + t17) - m(5) * pkin(3) - mrSges(5,1) * t29 + mrSges(5,2) * t28;
t50 = m(5) + m(6) + m(7);
t73 = -m(4) - t50;
t74 = -m(3) + t73;
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t68 = -g(1) * t32 + g(2) * t34;
t72 = t75 * t31 + mrSges(2,2) - mrSges(3,3) - t77;
t18 = sin(t27);
t60 = pkin(5) * t18;
t61 = pkin(4) * t28;
t70 = m(6) * t61 + m(7) * (t60 + t61) + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t28 * mrSges(5,1) + t29 * mrSges(5,2);
t69 = -m(7) * pkin(5) - mrSges(6,1);
t20 = qJ(6) + t27;
t15 = sin(t20);
t16 = cos(t20);
t67 = mrSges(6,1) * t19 + mrSges(7,1) * t16 - mrSges(6,2) * t18 - mrSges(7,2) * t15 - t75;
t55 = t32 * t15;
t5 = t16 * t34 - t31 * t55;
t54 = t32 * t16;
t6 = t15 * t34 + t31 * t54;
t63 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t56 = t31 * t34;
t7 = t15 * t56 + t54;
t8 = t16 * t56 - t55;
t62 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t57 = g(3) * t33;
t53 = t32 * t18;
t52 = t32 * t19;
t51 = t34 * pkin(1) + t32 * qJ(2);
t41 = -mrSges(7,1) * t15 - mrSges(7,2) * t16;
t11 = t18 * t56 + t52;
t9 = t19 * t34 - t31 * t53;
t12 = t19 * t56 - t53;
t10 = t18 * t34 + t31 * t52;
t1 = [(-m(3) * t51 - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t73 * (t34 * pkin(7) + t51) - t70 * t34 + t72 * t32) * g(2) + (-t12 * mrSges(6,1) - t8 * mrSges(7,1) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + (m(3) * pkin(1) + t73 * (-pkin(1) - pkin(7)) + t70) * t32 + (t74 * qJ(2) + t72) * t34) * g(1), -t68 * t74 (t67 * t31 + t77) * g(3) + t68 * ((mrSges(4,1) + t67) * t33 - t78 * t31) (-t31 * g(3) - t68 * t33) * t50 (m(7) * t60 + mrSges(6,1) * t18 + mrSges(6,2) * t19 - t41) * t57 + (-t12 * mrSges(6,2) + t69 * t11 - t62) * g(2) + (t10 * mrSges(6,2) + t69 * t9 - t63) * g(1), -g(1) * t63 - g(2) * t62 - t41 * t57];
taug  = t1(:);

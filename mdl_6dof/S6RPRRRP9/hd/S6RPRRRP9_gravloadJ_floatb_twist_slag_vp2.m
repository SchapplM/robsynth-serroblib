% Calculate Gravitation load on the joints for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:10
% EndTime: 2019-03-09 06:26:11
% DurationCPUTime: 0.77s
% Computational Cost: add. (333->88), mult. (494->102), div. (0->0), fcn. (446->8), ass. (0->49)
t37 = -pkin(9) - pkin(8);
t91 = mrSges(4,2) - m(5) * pkin(8) + m(6) * t37 + m(7) * (-qJ(6) + t37) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t90 = t32 * mrSges(4,1) + t91 * t35;
t30 = qJ(4) + qJ(5);
t23 = cos(t30);
t34 = cos(qJ(4));
t26 = t34 * pkin(4);
t19 = pkin(5) * t23 + t26;
t89 = -m(5) * pkin(3) - m(6) * (t26 + pkin(3)) - m(7) * (pkin(3) + t19);
t87 = mrSges(6,1) + mrSges(7,1);
t82 = mrSges(6,2) + mrSges(7,2);
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t74 = -g(1) * t33 + g(2) * t36;
t86 = m(6) * pkin(4);
t31 = sin(qJ(4));
t62 = pkin(4) * t31;
t85 = m(6) * t62;
t22 = sin(t30);
t18 = pkin(5) * t22 + t62;
t84 = m(7) * t18;
t81 = -m(4) - m(5) - m(6) - m(7);
t80 = m(3) - t81;
t79 = mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t84;
t78 = t89 * t32 + mrSges(2,2) - mrSges(3,3) - t90;
t76 = -mrSges(5,1) - t86;
t75 = mrSges(6,1) * t22 + t82 * t23;
t55 = t32 * t36;
t11 = t22 * t55 + t23 * t33;
t12 = -t22 * t33 + t23 * t55;
t73 = -t87 * t11 - t82 * t12;
t56 = t32 * t33;
t10 = t22 * t36 + t23 * t56;
t9 = -t22 * t56 + t23 * t36;
t72 = t82 * t10 - t87 * t9;
t71 = t34 * mrSges(5,1) - t31 * mrSges(5,2) - t82 * t22 + t87 * t23 - t89;
t68 = m(7) * pkin(5);
t59 = g(3) * t35;
t57 = t31 * t36;
t54 = t33 * t34;
t53 = t34 * t36;
t52 = t36 * pkin(1) + t33 * qJ(2);
t15 = t31 * t55 + t54;
t13 = -t31 * t56 + t53;
t16 = -t31 * t33 + t32 * t53;
t14 = t32 * t54 + t57;
t1 = [(-t57 * t86 - m(3) * t52 - t14 * mrSges(5,1) - t13 * mrSges(5,2) - t82 * t9 - t87 * t10 + t81 * (t36 * pkin(7) + t52) - t79 * t36 + t78 * t33) * g(2) + (-t16 * mrSges(5,1) + t15 * mrSges(5,2) - t87 * t12 + t82 * t11 + (m(3) * pkin(1) + t85 + t81 * (-pkin(1) - pkin(7)) + t79) * t33 + (-t80 * qJ(2) + t78) * t36) * g(1), t74 * t80 (t71 * t32 + t90) * g(3) - t74 * ((-mrSges(4,1) - t71) * t35 + t91 * t32) (mrSges(5,1) * t31 + mrSges(7,1) * t22 + mrSges(5,2) * t34 + t75 + t84 + t85) * t59 + (-t16 * mrSges(5,2) - m(7) * (t18 * t55 + t19 * t33) + t76 * t15 + t73) * g(2) + (t14 * mrSges(5,2) - m(7) * (-t18 * t56 + t19 * t36) + t76 * t13 + t72) * g(1) (-(-mrSges(7,1) - t68) * t22 + t75) * t59 + (-t11 * t68 + t73) * g(2) + (-t9 * t68 + t72) * g(1) (-g(3) * t32 - t74 * t35) * m(7)];
taug  = t1(:);

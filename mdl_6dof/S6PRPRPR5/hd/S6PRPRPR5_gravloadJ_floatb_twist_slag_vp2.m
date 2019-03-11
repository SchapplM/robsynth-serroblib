% Calculate Gravitation load on the joints for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:00
% EndTime: 2019-03-08 19:43:01
% DurationCPUTime: 0.73s
% Computational Cost: add. (458->87), mult. (786->123), div. (0->0), fcn. (881->12), ass. (0->50)
t84 = m(6) + m(7);
t91 = -mrSges(5,1) + mrSges(6,2);
t90 = mrSges(5,2) - mrSges(6,3);
t39 = sin(qJ(6));
t41 = cos(qJ(6));
t73 = -m(4) * qJ(3) - m(7) * pkin(5) - t41 * mrSges(7,1) + t39 * mrSges(7,2) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t33 = pkin(11) + qJ(4);
t31 = sin(t33);
t32 = cos(t33);
t37 = cos(pkin(11));
t89 = -t90 * t31 - t91 * t32 + mrSges(3,1) + m(4) * pkin(2) + t37 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2);
t88 = t39 * mrSges(7,1) + t41 * mrSges(7,2);
t78 = -m(7) * pkin(9) - mrSges(7,3);
t85 = pkin(4) * t84 - t78 - t91;
t80 = m(5) + t84;
t76 = -t84 * qJ(5) - t88 + t90;
t74 = t32 * mrSges(7,3) + t88 * t31 + t89;
t35 = sin(pkin(10));
t40 = sin(qJ(2));
t42 = cos(qJ(2));
t56 = cos(pkin(10));
t57 = cos(pkin(6));
t46 = t57 * t56;
t18 = t35 * t40 - t42 * t46;
t19 = t35 * t42 + t40 * t46;
t30 = pkin(3) * t37 + pkin(2);
t38 = -pkin(8) - qJ(3);
t71 = -t18 * t30 - t19 * t38;
t70 = t18 * t32;
t51 = t35 * t57;
t20 = t40 * t56 + t42 * t51;
t69 = t20 * t32;
t36 = sin(pkin(6));
t64 = t35 * t36;
t63 = t36 * t40;
t62 = t36 * t42;
t61 = t39 * t42;
t60 = t41 * t42;
t21 = -t40 * t51 + t42 * t56;
t59 = -t20 * t30 - t21 * t38;
t58 = qJ(5) * t31;
t54 = m(4) + t80;
t53 = -pkin(4) * t70 - t18 * t58 + t71;
t50 = t36 * t56;
t49 = -pkin(4) * t69 - t20 * t58 + t59;
t23 = t30 * t62;
t14 = t31 * t63 - t32 * t57;
t5 = t21 * t31 - t32 * t64;
t3 = t19 * t31 + t32 * t50;
t1 = [(-m(2) - m(3) - t54) * g(3) (-m(5) * t71 - m(6) * t53 - m(7) * (-pkin(9) * t70 + t53) + t73 * t19 + t74 * t18) * g(2) + (-m(5) * t59 - m(6) * t49 - m(7) * (-pkin(9) * t69 + t49) + t73 * t21 + t74 * t20) * g(1) + (-m(5) * t23 - t84 * (t23 + (pkin(4) * t32 + t58) * t62) + ((-mrSges(7,1) * t61 - mrSges(7,2) * t60) * t31 + (t80 * t38 + t73) * t40 + (t78 * t32 - t89) * t42) * t36) * g(3) (-g(1) * t20 - g(2) * t18 + g(3) * t62) * t54 (t76 * (t31 * t57 + t32 * t63) + t85 * t14) * g(3) + (t76 * (t19 * t32 - t31 * t50) + t85 * t3) * g(2) + (t76 * (t21 * t32 + t31 * t64) + t85 * t5) * g(1), t84 * (-g(1) * t5 - g(2) * t3 - g(3) * t14) -g(1) * ((-t20 * t39 + t41 * t5) * mrSges(7,1) + (-t20 * t41 - t39 * t5) * mrSges(7,2)) - g(2) * ((-t18 * t39 + t3 * t41) * mrSges(7,1) + (-t18 * t41 - t3 * t39) * mrSges(7,2)) - g(3) * ((t14 * t41 + t36 * t61) * mrSges(7,1) + (-t14 * t39 + t36 * t60) * mrSges(7,2))];
taug  = t1(:);

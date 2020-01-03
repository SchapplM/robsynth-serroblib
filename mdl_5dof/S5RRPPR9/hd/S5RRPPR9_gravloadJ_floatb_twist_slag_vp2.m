% Calculate Gravitation load on the joints for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:41
% EndTime: 2019-12-31 19:40:43
% DurationCPUTime: 0.50s
% Computational Cost: add. (147->71), mult. (309->79), div. (0->0), fcn. (262->6), ass. (0->38)
t69 = -mrSges(3,1) - mrSges(4,1);
t68 = mrSges(5,2) - mrSges(6,3);
t67 = mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t53 = -pkin(2) - pkin(3);
t66 = -m(5) * t53 - m(6) * (-pkin(7) + t53) - t68;
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t59 = g(1) * t21 + g(2) * t18;
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t63 = t67 * t17 - t69 * t20;
t61 = m(5) + m(6);
t58 = m(4) + t61;
t56 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t55 = t68 * t20 - t63;
t50 = g(3) * t20;
t13 = t20 * pkin(2);
t49 = t20 * pkin(7);
t16 = sin(qJ(5));
t47 = t18 * t16;
t19 = cos(qJ(5));
t46 = t18 * t19;
t44 = t21 * t16;
t43 = t21 * t19;
t42 = t21 * t20;
t10 = t17 * qJ(3);
t41 = t13 + t10;
t40 = t21 * pkin(1) + t18 * pkin(6);
t37 = t20 * pkin(3) + t41;
t36 = pkin(2) * t42 + t21 * t10 + t40;
t35 = -pkin(1) - t10;
t25 = m(6) * pkin(4) + t19 * mrSges(6,1) - t16 * mrSges(6,2);
t14 = t21 * pkin(6);
t4 = t17 * t43 - t47;
t3 = -t17 * t44 - t46;
t2 = -t17 * t46 - t44;
t1 = t17 * t47 - t43;
t5 = [(-m(3) * t40 - m(4) * t36 - t4 * mrSges(6,1) - t3 * mrSges(6,2) - t61 * (pkin(3) * t42 - t18 * qJ(4) + t36) + t56 * t18 + (-mrSges(2,1) - m(6) * (t17 * pkin(4) + t49) + t55) * t21) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) - t61 * (-t21 * qJ(4) + t14) + (-m(3) - m(4)) * t14 + t56 * t21 + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t35 - t13) - m(5) * t35 - m(6) * (-pkin(1) + (-pkin(4) - qJ(3)) * t17) + t66 * t20 + t63) * t18) * g(1), (-m(4) * t41 - m(5) * t37 - m(6) * (t37 + t49) - t25 * t17 + t55) * g(3) + ((m(4) * pkin(2) + t66 - t69) * t17 + (-qJ(3) * t58 - t25 - t67) * t20) * t59, (-t59 * t17 + t50) * t58, t61 * (g(1) * t18 - g(2) * t21), -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - (mrSges(6,1) * t16 + mrSges(6,2) * t19) * t50];
taug = t5(:);

% Calculate Gravitation load on the joints for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:27
% DurationCPUTime: 0.34s
% Computational Cost: add. (267->56), mult. (194->57), div. (0->0), fcn. (147->8), ass. (0->30)
t63 = mrSges(5,1) + mrSges(6,1);
t62 = -mrSges(5,2) + mrSges(6,3);
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t61 = t62 * t24 + t63 * t26;
t60 = -mrSges(4,1) - t61;
t59 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t33 = pkin(4) * t26 + qJ(5) * t24;
t23 = qJ(1) + qJ(2);
t19 = pkin(8) + t23;
t16 = sin(t19);
t17 = cos(t19);
t57 = g(1) * t17 + g(2) * t16;
t20 = sin(t23);
t21 = cos(t23);
t56 = mrSges(3,1) * t20 + mrSges(3,2) * t21 + t59 * t17 + (-m(6) * (-pkin(3) - t33) - t60) * t16;
t55 = -t21 * mrSges(3,1) + t20 * mrSges(3,2) + t59 * t16 + t60 * t17;
t25 = sin(qJ(1));
t54 = pkin(1) * t25;
t53 = pkin(2) * t20;
t18 = pkin(2) * t21;
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t43 = t17 * pkin(3) + t16 * pkin(7) + t18;
t13 = t17 * pkin(7);
t42 = t13 - t53;
t39 = t33 * t17 + t43;
t38 = -t53 - t54;
t30 = -pkin(3) * t16 + t42;
t1 = [(-mrSges(2,1) * t27 + t25 * mrSges(2,2) - m(3) * t22 - m(4) * (t18 + t22) - m(5) * (t22 + t43) - m(6) * (t22 + t39) + t55) * g(2) + (t25 * mrSges(2,1) + mrSges(2,2) * t27 + m(3) * t54 - m(4) * t38 - m(5) * (t30 - t54) - m(6) * (t13 + t38) + t56) * g(1), (-m(4) * t18 - m(5) * t43 - m(6) * t39 + t55) * g(2) + (m(4) * t53 - m(5) * t30 - m(6) * t42 + t56) * g(1), (-m(4) - m(5) - m(6)) * g(3), (-m(6) * t33 - t61) * g(3) + ((-m(6) * qJ(5) - t62) * t26 + (m(6) * pkin(4) + t63) * t24) * t57, (g(3) * t26 - t57 * t24) * m(6)];
taug = t1(:);

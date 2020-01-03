% Calculate Gravitation load on the joints for
% S5RRRRP3
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:07
% DurationCPUTime: 0.32s
% Computational Cost: add. (348->61), mult. (241->62), div. (0->0), fcn. (187->8), ass. (0->34)
t67 = mrSges(5,1) + mrSges(6,1);
t66 = -mrSges(5,2) + mrSges(6,3);
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t65 = t66 * t24 + t67 * t26;
t64 = -mrSges(4,1) - t65;
t63 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t33 = t26 * pkin(4) + t24 * qJ(5);
t23 = qJ(1) + qJ(2);
t21 = qJ(3) + t23;
t17 = sin(t21);
t18 = cos(t21);
t61 = g(1) * t18 + g(2) * t17;
t60 = t63 * t18 + (-m(6) * (-pkin(3) - t33) - t64) * t17;
t59 = t63 * t17 + t64 * t18;
t19 = sin(t23);
t20 = cos(t23);
t58 = t19 * mrSges(3,1) + t20 * mrSges(3,2) + t60;
t57 = -t20 * mrSges(3,1) + t19 * mrSges(3,2) + t59;
t56 = pkin(2) * t19;
t16 = pkin(2) * t20;
t25 = sin(qJ(1));
t53 = t25 * pkin(1);
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t47 = t18 * pkin(3) + t17 * pkin(8);
t44 = t16 + t47;
t13 = t18 * pkin(8);
t43 = -t17 * pkin(3) + t13;
t42 = t18 * t33 + t47;
t39 = t16 + t42;
t38 = -t53 - t56;
t30 = t43 - t56;
t1 = [(-t27 * mrSges(2,1) + t25 * mrSges(2,2) - m(3) * t22 - m(4) * (t16 + t22) - m(5) * (t22 + t44) - m(6) * (t22 + t39) + t57) * g(2) + (t25 * mrSges(2,1) + t27 * mrSges(2,2) + m(3) * t53 - m(4) * t38 - m(5) * (t30 - t53) - m(6) * (t13 + t38) + t58) * g(1), (-m(4) * t16 - m(5) * t44 - m(6) * t39 + t57) * g(2) + (m(4) * t56 - m(5) * t30 - m(6) * (t13 - t56) + t58) * g(1), (-m(5) * t47 - m(6) * t42 + t59) * g(2) + (-m(5) * t43 - m(6) * t13 + t60) * g(1), (-m(6) * t33 - t65) * g(3) + ((-m(6) * qJ(5) - t66) * t26 + (m(6) * pkin(4) + t67) * t24) * t61, (g(3) * t26 - t61 * t24) * m(6)];
taug = t1(:);

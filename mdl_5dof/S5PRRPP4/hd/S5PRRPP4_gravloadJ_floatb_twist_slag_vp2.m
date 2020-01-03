% Calculate Gravitation load on the joints for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:46
% EndTime: 2019-12-31 17:40:48
% DurationCPUTime: 0.34s
% Computational Cost: add. (161->43), mult. (177->42), div. (0->0), fcn. (135->4), ass. (0->20)
t44 = -mrSges(4,1) - mrSges(5,1);
t43 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t28 = m(5) + m(6);
t38 = -m(4) - t28;
t12 = pkin(7) + qJ(2);
t8 = sin(t12);
t9 = cos(t12);
t35 = g(1) * t9 + g(2) * t8;
t41 = g(1) * t8;
t13 = sin(qJ(3));
t10 = t13 * qJ(4);
t14 = cos(qJ(3));
t11 = t14 * pkin(3);
t25 = t11 + t10;
t34 = t43 * t13 + t44 * t14;
t32 = -mrSges(3,1) + t34;
t31 = m(6) * qJ(5) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t22 = -pkin(2) - t10;
t20 = m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1);
t1 = [(-m(2) - m(3) + t38) * g(3), (t31 * t8 + t38 * (t9 * pkin(2) + t8 * pkin(6))) * g(2) + (m(4) * pkin(2) - m(5) * (t22 - t11) - m(6) * t22 - t20 * t14 - t32) * t41 + ((-t28 * t25 - (m(6) * pkin(4) + mrSges(6,1)) * t14 + t32) * g(2) + (t38 * pkin(6) + t31) * g(1)) * t9, (-m(5) * t25 - m(6) * (pkin(4) * t14 + t25) - t14 * mrSges(6,1) + t34) * g(3) + ((m(5) * pkin(3) - t20 - t44) * t13 + (-qJ(4) * t28 + t43) * t14) * t35, (g(3) * t14 - t13 * t35) * t28, (-g(2) * t9 + t41) * m(6)];
taug = t1(:);

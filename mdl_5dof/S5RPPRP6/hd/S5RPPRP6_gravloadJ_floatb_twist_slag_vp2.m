% Calculate Gravitation load on the joints for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:54:59
% DurationCPUTime: 0.25s
% Computational Cost: add. (123->38), mult. (169->41), div. (0->0), fcn. (125->6), ass. (0->20)
t39 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t38 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t9 = pkin(7) + qJ(4);
t4 = sin(t9);
t5 = cos(t9);
t37 = -t38 * t5 + t39 * t4;
t13 = sin(qJ(1));
t14 = cos(qJ(1));
t35 = -g(1) * t13 + g(2) * t14;
t36 = -m(5) - m(6);
t33 = mrSges(2,1) + mrSges(6,2) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t10 = sin(pkin(7));
t32 = mrSges(2,2) - mrSges(3,3) - t10 * mrSges(4,1) - cos(pkin(7)) * mrSges(4,2) - t37;
t31 = t14 * pkin(1) + t13 * qJ(2);
t30 = pkin(3) * t10;
t27 = -m(4) + t36;
t7 = t14 * qJ(2);
t25 = -t13 * pkin(1) + t7;
t12 = -pkin(6) - qJ(3);
t1 = [((-m(3) - m(4)) * t31 + t36 * (-t12 * t14 + t13 * t30 + t31) + (-m(4) * qJ(3) - t33) * t14 + t32 * t13) * g(2) + (-m(3) * t25 - m(4) * t7 + t36 * (t13 * t12 + t14 * t30 + t25) + t32 * t14 + (-m(4) * (-pkin(1) - qJ(3)) + t33) * t13) * g(1), t35 * (m(3) - t27), (g(1) * t14 + g(2) * t13) * t27, t37 * g(3) + (t38 * t4 + t39 * t5) * t35, (-g(3) * t4 - t35 * t5) * m(6)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (116->35), mult. (125->35), div. (0->0), fcn. (93->8), ass. (0->19)
t11 = pkin(7) + qJ(3);
t8 = qJ(4) + t11;
t3 = sin(t8);
t4 = cos(t8);
t25 = t4 * mrSges(5,1) - t3 * mrSges(5,2);
t6 = sin(t11);
t35 = -t6 * mrSges(4,2) + t25;
t15 = sin(qJ(1));
t16 = cos(qJ(1));
t34 = g(1) * t16 + g(2) * t15;
t13 = cos(pkin(7));
t5 = t13 * pkin(2) + pkin(1);
t7 = cos(t11);
t33 = mrSges(2,1) + m(3) * pkin(1) + t13 * mrSges(3,1) - sin(pkin(7)) * mrSges(3,2) + m(5) * (pkin(3) * t7 + t5) + m(4) * t5 + t7 * mrSges(4,1) + t35;
t14 = -pkin(5) - qJ(2);
t32 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) + m(5) * (-pkin(6) + t14) - mrSges(5,3) + m(4) * t14 - mrSges(4,3);
t26 = m(5) * pkin(3) + mrSges(4,1);
t21 = mrSges(5,1) * t3 + mrSges(5,2) * t4;
t1 = [(t32 * t15 - t33 * t16) * g(2) + (t33 * t15 + t32 * t16) * g(1), (-g(1) * t15 + g(2) * t16) * (m(3) + m(4) + m(5)), (-t26 * t7 - t35) * g(3) + t34 * (mrSges(4,2) * t7 + t26 * t6 + t21), -g(3) * t25 + t34 * t21];
taug = t1(:);

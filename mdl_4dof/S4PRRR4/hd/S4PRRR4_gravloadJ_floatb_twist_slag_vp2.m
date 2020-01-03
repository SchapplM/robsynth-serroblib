% Calculate Gravitation load on the joints for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:30
% DurationCPUTime: 0.16s
% Computational Cost: add. (103->28), mult. (95->28), div. (0->0), fcn. (69->6), ass. (0->16)
t8 = qJ(3) + qJ(4);
t5 = sin(t8);
t6 = cos(t8);
t17 = t6 * mrSges(5,1) - t5 * mrSges(5,2);
t9 = sin(qJ(3));
t27 = -t9 * mrSges(4,2) + t17;
t7 = pkin(7) + qJ(2);
t3 = sin(t7);
t4 = cos(t7);
t26 = g(1) * t4 + g(2) * t3;
t10 = cos(qJ(3));
t25 = m(4) * pkin(2) + t10 * mrSges(4,1) + mrSges(3,1) + m(5) * (pkin(3) * t10 + pkin(2)) + t27;
t24 = mrSges(3,2) + m(5) * (-pkin(6) - pkin(5)) - mrSges(5,3) - m(4) * pkin(5) - mrSges(4,3);
t18 = m(5) * pkin(3) + mrSges(4,1);
t15 = mrSges(5,1) * t5 + mrSges(5,2) * t6;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), (t24 * t3 - t25 * t4) * g(2) + (t24 * t4 + t25 * t3) * g(1), (-t18 * t10 - t27) * g(3) + t26 * (mrSges(4,2) * t10 + t18 * t9 + t15), -g(3) * t17 + t26 * t15];
taug = t1(:);

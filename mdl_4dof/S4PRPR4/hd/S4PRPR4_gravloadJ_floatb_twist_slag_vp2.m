% Calculate Gravitation load on the joints for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (66->20), mult. (63->23), div. (0->0), fcn. (42->4), ass. (0->11)
t18 = m(4) + m(5);
t6 = pkin(6) + qJ(2);
t4 = sin(t6);
t5 = cos(t6);
t17 = -g(1) * t4 + g(2) * t5;
t7 = sin(qJ(4));
t8 = cos(qJ(4));
t9 = t7 * mrSges(5,1) + t8 * mrSges(5,2);
t16 = mrSges(3,2) - mrSges(4,3) - t9;
t15 = mrSges(3,1) + mrSges(5,3) - mrSges(4,2);
t1 = [(-m(2) - m(3) - t18) * g(3), (-t18 * (t5 * pkin(2) + t4 * qJ(3)) + (-m(5) * pkin(5) - t15) * t5 + t16 * t4) * g(2) + ((m(4) * pkin(2) - m(5) * (-pkin(2) - pkin(5)) + t15) * t4 + (-t18 * qJ(3) + t16) * t5) * g(1), t18 * t17, g(3) * t9 + t17 * (mrSges(5,1) * t8 - mrSges(5,2) * t7)];
taug = t1(:);

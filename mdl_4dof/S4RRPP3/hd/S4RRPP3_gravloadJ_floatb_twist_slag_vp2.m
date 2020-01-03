% Calculate Gravitation load on the joints for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:23
% DurationCPUTime: 0.25s
% Computational Cost: add. (112->39), mult. (153->41), div. (0->0), fcn. (117->6), ass. (0->19)
t30 = m(4) + m(5);
t34 = -mrSges(4,1) - mrSges(5,1);
t33 = mrSges(4,2) - mrSges(5,3);
t11 = cos(qJ(1));
t9 = sin(qJ(1));
t32 = g(1) * t11 + g(2) * t9;
t10 = cos(qJ(2));
t6 = qJ(2) + pkin(6);
t3 = sin(t6);
t4 = cos(t6);
t8 = sin(qJ(2));
t31 = -t10 * mrSges(3,1) + t8 * mrSges(3,2) + t33 * t3 + t34 * t4;
t27 = m(3) * pkin(1) + mrSges(2,1) - t31;
t26 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t5 = t10 * pkin(2);
t17 = pkin(3) * t4 + qJ(4) * t3;
t7 = -qJ(3) - pkin(5);
t2 = t5 + pkin(1);
t1 = [(-t30 * (t11 * t2 - t9 * t7) + t26 * t9 + (-m(5) * t17 - t27) * t11) * g(2) + ((m(4) * t2 - m(5) * (-t17 - t2) + t27) * t9 + (t30 * t7 + t26) * t11) * g(1), (-m(4) * t5 - m(5) * (t17 + t5) + t31) * g(3) + t32 * (mrSges(3,2) * t10 + (-m(5) * qJ(4) + t33) * t4 + (m(5) * pkin(3) - t34) * t3 + (t30 * pkin(2) + mrSges(3,1)) * t8), t30 * (-g(1) * t9 + g(2) * t11), (g(3) * t4 - t32 * t3) * m(5)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (90->27), mult. (99->27), div. (0->0), fcn. (71->6), ass. (0->16)
t14 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t15 = mrSges(4,2) + mrSges(5,2);
t7 = sin(qJ(3));
t9 = cos(qJ(3));
t24 = -t14 * t9 + t15 * t7;
t23 = m(4) + m(5);
t19 = m(3) + t23;
t20 = pkin(1) * t19 + mrSges(2,1);
t18 = t23 * pkin(2) + mrSges(3,1) - t24;
t17 = -m(4) * pkin(5) + m(5) * (-qJ(4) - pkin(5)) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t10 = cos(qJ(1));
t8 = sin(qJ(1));
t5 = qJ(1) + pkin(6);
t3 = cos(t5);
t2 = sin(t5);
t1 = [(t8 * mrSges(2,2) - t20 * t10 + t17 * t2 - t18 * t3) * g(2) + (mrSges(2,2) * t10 + t17 * t3 + t18 * t2 + t20 * t8) * g(1), -t19 * g(3), t24 * g(3) + (g(1) * t3 + g(2) * t2) * (t14 * t7 + t15 * t9), (-g(1) * t2 + g(2) * t3) * m(5)];
taug = t1(:);

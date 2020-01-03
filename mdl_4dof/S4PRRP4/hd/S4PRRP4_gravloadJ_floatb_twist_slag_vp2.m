% Calculate Gravitation load on the joints for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:45
% EndTime: 2019-12-31 16:27:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (91->27), mult. (101->29), div. (0->0), fcn. (75->4), ass. (0->14)
t24 = mrSges(4,1) + mrSges(5,1);
t23 = -mrSges(4,2) + mrSges(5,3);
t7 = sin(qJ(3));
t8 = cos(qJ(3));
t22 = t23 * t7 + t24 * t8;
t21 = -m(4) - m(5);
t11 = pkin(3) * t8 + t7 * qJ(4);
t20 = -m(5) * t11 - t22;
t6 = pkin(6) + qJ(2);
t4 = sin(t6);
t5 = cos(t6);
t19 = g(1) * t5 + g(2) * t4;
t18 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t1 = [(-m(2) - m(3) + t21) * g(3), (t21 * (t5 * pkin(2) + t4 * pkin(5)) + (-mrSges(3,1) + t20) * t5 + t18 * t4) * g(2) + ((mrSges(3,1) + m(4) * pkin(2) - m(5) * (-pkin(2) - t11) + t22) * t4 + (t21 * pkin(5) + t18) * t5) * g(1), t20 * g(3) + ((-m(5) * qJ(4) - t23) * t8 + (m(5) * pkin(3) + t24) * t7) * t19, (g(3) * t8 - t19 * t7) * m(5)];
taug = t1(:);

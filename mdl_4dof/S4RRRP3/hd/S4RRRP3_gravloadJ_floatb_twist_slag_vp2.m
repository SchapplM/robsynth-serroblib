% Calculate Gravitation load on the joints for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:02
% DurationCPUTime: 0.22s
% Computational Cost: add. (153->41), mult. (161->46), div. (0->0), fcn. (125->6), ass. (0->23)
t53 = mrSges(4,1) + mrSges(5,1);
t52 = -mrSges(4,2) + mrSges(5,3);
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t51 = t52 * t19 + t53 * t21;
t50 = -mrSges(3,1) - t51;
t49 = -mrSges(5,2) - mrSges(4,3) + mrSges(3,2);
t27 = pkin(3) * t21 + qJ(4) * t19;
t18 = qJ(1) + qJ(2);
t15 = sin(t18);
t16 = cos(t18);
t47 = g(1) * t16 + g(2) * t15;
t46 = t49 * t16 + (-m(5) * (-pkin(2) - t27) - t50) * t15;
t45 = t49 * t15 + t50 * t16;
t20 = sin(qJ(1));
t44 = pkin(1) * t20;
t22 = cos(qJ(1));
t17 = t22 * pkin(1);
t36 = t16 * pkin(2) + t15 * pkin(6);
t13 = t16 * pkin(6);
t33 = -pkin(2) * t15 + t13;
t32 = t27 * t16 + t36;
t1 = [(-mrSges(2,1) * t22 + t20 * mrSges(2,2) - m(3) * t17 - m(4) * (t17 + t36) - m(5) * (t17 + t32) + t45) * g(2) + (t20 * mrSges(2,1) + mrSges(2,2) * t22 + m(3) * t44 - m(4) * (t33 - t44) - m(5) * (t13 - t44) + t46) * g(1), (-m(4) * t36 - m(5) * t32 + t45) * g(2) + (-m(4) * t33 - m(5) * t13 + t46) * g(1), (-m(5) * t27 - t51) * g(3) + ((-m(5) * qJ(4) - t52) * t21 + (m(5) * pkin(3) + t53) * t19) * t47, (g(3) * t21 - t47 * t19) * m(5)];
taug = t1(:);

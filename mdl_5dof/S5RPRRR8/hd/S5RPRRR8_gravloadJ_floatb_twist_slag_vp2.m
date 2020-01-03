% Calculate Gravitation load on the joints for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:42
% EndTime: 2019-12-31 19:05:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (192->44), mult. (313->50), div. (0->0), fcn. (327->8), ass. (0->26)
t54 = -mrSges(4,2) - m(6) * (-pkin(8) - pkin(7)) + mrSges(6,3) + m(5) * pkin(7) + mrSges(5,3);
t16 = cos(qJ(4));
t15 = sin(qJ(4));
t14 = qJ(4) + qJ(5);
t8 = sin(t14);
t9 = cos(t14);
t27 = -t9 * mrSges(6,1) + t8 * mrSges(6,2);
t50 = -t15 * mrSges(5,2) - t27;
t53 = mrSges(4,1) + m(5) * pkin(3) + t16 * mrSges(5,1) + m(6) * (t16 * pkin(4) + pkin(3)) + t50;
t33 = sin(qJ(3));
t34 = sin(qJ(1));
t35 = cos(qJ(3));
t36 = cos(qJ(1));
t1 = -t34 * t33 - t36 * t35;
t2 = t36 * t33 - t34 * t35;
t52 = t54 * t1 + t53 * t2;
t51 = -t53 * t1 + t54 * t2;
t47 = mrSges(2,1) + mrSges(3,1);
t46 = mrSges(2,2) - mrSges(3,3);
t43 = g(1) * t1 + g(2) * t2;
t42 = -m(4) - m(5) - m(6);
t31 = t36 * pkin(1) + t34 * qJ(2);
t28 = m(6) * pkin(4) + mrSges(5,1);
t25 = -t34 * pkin(1) + t36 * qJ(2);
t22 = mrSges(6,1) * t8 + mrSges(6,2) * t9;
t3 = [(-m(3) * t31 - t47 * t36 + t46 * t34 + t42 * (t36 * pkin(2) + t31) - t51) * g(2) + (-m(3) * t25 + t46 * t36 + t47 * t34 + t42 * (-t34 * pkin(2) + t25) - t52) * g(1), (-t34 * g(1) + t36 * g(2)) * (m(3) - t42), t52 * g(1) + t51 * g(2), (t28 * t16 + t50) * g(3) + t43 * (-mrSges(5,2) * t16 - t28 * t15 - t22), -g(3) * t27 - t22 * t43];
taug = t3(:);

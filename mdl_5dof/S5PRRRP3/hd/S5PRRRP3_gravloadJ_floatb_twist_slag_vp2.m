% Calculate Gravitation load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:31
% EndTime: 2019-12-05 16:43:31
% DurationCPUTime: 0.23s
% Computational Cost: add. (192->44), mult. (167->41), div. (0->0), fcn. (122->6), ass. (0->22)
t15 = sin(qJ(3));
t14 = qJ(3) + qJ(4);
t10 = cos(t14);
t29 = mrSges(5,2) + mrSges(6,2);
t30 = mrSges(5,1) + mrSges(6,1);
t9 = sin(t14);
t23 = -t30 * t10 + t29 * t9;
t40 = -t15 * mrSges(4,2) - t23;
t13 = pkin(8) + qJ(2);
t7 = sin(t13);
t8 = cos(t13);
t37 = g(1) * t8 + g(2) * t7;
t16 = cos(qJ(3));
t11 = t16 * pkin(3);
t5 = pkin(4) * t10;
t32 = t5 + t11;
t36 = m(4) * pkin(2) + t16 * mrSges(4,1) + mrSges(3,1) + m(6) * (pkin(2) + t32) + m(5) * (t11 + pkin(2)) + t40;
t17 = -pkin(7) - pkin(6);
t35 = mrSges(3,2) + m(6) * (-qJ(5) + t17) - mrSges(6,3) + m(5) * t17 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t27 = m(5) * pkin(3) + mrSges(4,1);
t24 = t29 * t10;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), (t35 * t7 - t36 * t8) * g(2) + (t35 * t8 + t36 * t7) * g(1), (-m(6) * t32 - t27 * t16 - t40) * g(3) + t37 * (-m(6) * (-pkin(3) * t15 - pkin(4) * t9) + mrSges(4,2) * t16 + t27 * t15 + t30 * t9 + t24), (-m(6) * t5 + t23) * g(3) + t37 * (t24 + (m(6) * pkin(4) + t30) * t9), (-g(1) * t7 + g(2) * t8) * m(6)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:59
% EndTime: 2019-12-05 15:03:00
% DurationCPUTime: 0.22s
% Computational Cost: add. (108->30), mult. (141->41), div. (0->0), fcn. (108->6), ass. (0->18)
t37 = mrSges(4,1) + mrSges(6,3) - mrSges(5,2);
t10 = sin(qJ(5));
t11 = cos(qJ(5));
t36 = -t10 * mrSges(6,1) - t11 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t29 = g(1) * t9 + g(2) * t8;
t24 = m(5) + m(6);
t7 = pkin(8) + qJ(3);
t6 = cos(t7);
t25 = g(3) * t6;
t22 = t11 * t8;
t21 = t11 * t9;
t20 = t8 * t10;
t19 = t9 * t10;
t17 = m(3) + m(4) + t24;
t5 = sin(t7);
t1 = [(-m(2) - t17) * g(3), (-g(1) * t8 + g(2) * t9) * t17, (-t24 * (t6 * pkin(3) + t5 * qJ(4)) + (-m(6) * pkin(6) - t37) * t6 + t36 * t5) * g(3) + ((-m(6) * (-pkin(3) - pkin(6)) + m(5) * pkin(3) + t37) * t5 + (-t24 * qJ(4) + t36) * t6) * t29, (-t29 * t5 + t25) * t24, -g(1) * ((t5 * t21 - t20) * mrSges(6,1) + (-t5 * t19 - t22) * mrSges(6,2)) - g(2) * ((t5 * t22 + t19) * mrSges(6,1) + (-t5 * t20 + t21) * mrSges(6,2)) - (-mrSges(6,1) * t11 + mrSges(6,2) * t10) * t25];
taug = t1(:);

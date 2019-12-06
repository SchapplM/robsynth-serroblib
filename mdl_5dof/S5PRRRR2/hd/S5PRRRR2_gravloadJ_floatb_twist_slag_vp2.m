% Calculate Gravitation load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:43
% EndTime: 2019-12-05 17:04:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (164->38), mult. (122->40), div. (0->0), fcn. (84->8), ass. (0->23)
t16 = qJ(2) + qJ(3);
t14 = qJ(4) + t16;
t10 = sin(t14);
t11 = cos(t14);
t26 = -m(6) * pkin(6) + mrSges(5,2);
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t34 = mrSges(6,1) * t19 - mrSges(6,2) * t17;
t33 = mrSges(5,1) + t34;
t36 = t33 * t10 + (-mrSges(6,3) + t26) * t11;
t31 = m(5) + m(6);
t13 = cos(t16);
t9 = pkin(3) * t13;
t20 = cos(qJ(2));
t30 = t20 * pkin(2) + t9;
t27 = m(4) + t31;
t24 = -t10 * mrSges(6,3) - t33 * t11;
t12 = sin(t16);
t22 = -t13 * mrSges(4,1) + t12 * mrSges(4,2) + t10 * mrSges(5,2) + t24;
t21 = mrSges(4,2) * t13 + (t31 * pkin(3) + mrSges(4,1)) * t12 + t36;
t18 = sin(qJ(2));
t6 = t10 * pkin(6);
t1 = [(-m(2) - m(3) - t27) * g(3), (t18 * mrSges(3,2) - m(5) * t30 - m(6) * (t6 + t30) + (-m(4) * pkin(2) - mrSges(3,1)) * t20 + t22) * g(2) + (mrSges(3,2) * t20 + (t27 * pkin(2) + mrSges(3,1)) * t18 + t21) * g(1), (-m(5) * t9 - m(6) * (t6 + t9) + t22) * g(2) + t21 * g(1), (t26 * t10 + t24) * g(2) + t36 * g(1), -g(3) * t34 + (g(1) * t11 + g(2) * t10) * (mrSges(6,1) * t17 + mrSges(6,2) * t19)];
taug = t1(:);

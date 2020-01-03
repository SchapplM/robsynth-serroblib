% Calculate Gravitation load on the joints for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (123->34), mult. (90->37), div. (0->0), fcn. (64->8), ass. (0->20)
t32 = mrSges(4,2) - mrSges(5,3);
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t31 = t19 * mrSges(5,1) - t17 * mrSges(5,2);
t30 = -mrSges(4,1) - t31;
t29 = m(4) + m(5);
t16 = qJ(1) + pkin(7);
t14 = qJ(3) + t16;
t10 = sin(t14);
t11 = cos(t14);
t28 = t11 * pkin(3) + t10 * pkin(6);
t13 = cos(t16);
t20 = cos(qJ(1));
t27 = t20 * pkin(1) + pkin(2) * t13;
t24 = m(3) + t29;
t22 = t32 * t10 + t30 * t11;
t21 = (m(5) * pkin(3) - t30) * t10 + (-m(5) * pkin(6) + t32) * t11;
t18 = sin(qJ(1));
t12 = sin(t16);
t1 = [(t18 * mrSges(2,2) - t13 * mrSges(3,1) + t12 * mrSges(3,2) - m(4) * t27 - m(5) * (t27 + t28) + (-m(3) * pkin(1) - mrSges(2,1)) * t20 + t22) * g(2) + (mrSges(2,2) * t20 + mrSges(3,2) * t13 + (t29 * pkin(2) + mrSges(3,1)) * t12 + (t24 * pkin(1) + mrSges(2,1)) * t18 + t21) * g(1), -t24 * g(3), (-m(5) * t28 + t22) * g(2) + t21 * g(1), -g(3) * t31 + (g(1) * t11 + g(2) * t10) * (mrSges(5,1) * t17 + mrSges(5,2) * t19)];
taug = t1(:);

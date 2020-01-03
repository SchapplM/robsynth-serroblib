% Calculate Gravitation load on the joints for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:32
% EndTime: 2019-12-31 16:51:33
% DurationCPUTime: 0.20s
% Computational Cost: add. (92->38), mult. (169->47), div. (0->0), fcn. (172->6), ass. (0->20)
t10 = sin(qJ(4));
t11 = cos(qJ(4));
t16 = -mrSges(5,1) * t11 + mrSges(5,2) * t10;
t31 = m(5) * pkin(3) + mrSges(4,1) - t16;
t30 = mrSges(2,1) + mrSges(3,1);
t29 = mrSges(2,2) - mrSges(3,3);
t12 = cos(qJ(1));
t22 = sin(qJ(3));
t23 = sin(qJ(1));
t24 = cos(qJ(3));
t1 = -t12 * t24 - t23 * t22;
t2 = t12 * t22 - t23 * t24;
t27 = t2 * mrSges(4,2) + t31 * t1;
t26 = t1 * mrSges(4,2) - t31 * t2;
t25 = t12 * pkin(1) + t23 * qJ(2);
t21 = t12 * pkin(2) + t25;
t20 = -m(5) * pkin(6) - mrSges(5,3);
t19 = -t23 * pkin(1) + t12 * qJ(2);
t13 = -t23 * pkin(2) + t19;
t3 = [(-m(3) * t25 - m(4) * t21 - m(5) * (pkin(6) * t2 + t21) - t2 * mrSges(5,3) + t29 * t23 - t30 * t12 + t27) * g(2) + (-m(3) * t19 - m(4) * t13 - m(5) * (t1 * pkin(6) + t13) - t1 * mrSges(5,3) + t30 * t23 + t29 * t12 + t26) * g(1), (-t23 * g(1) + g(2) * t12) * (m(3) + m(4) + m(5)), (-t20 * t2 - t27) * g(2) + (-t20 * t1 - t26) * g(1), -g(3) * t16 + (-g(1) * t1 - g(2) * t2) * (mrSges(5,1) * t10 + mrSges(5,2) * t11)];
taug = t3(:);

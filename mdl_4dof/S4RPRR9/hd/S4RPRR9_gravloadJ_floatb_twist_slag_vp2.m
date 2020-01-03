% Calculate Gravitation load on the joints for
% S4RPRR9
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:08
% EndTime: 2019-12-31 16:56:09
% DurationCPUTime: 0.23s
% Computational Cost: add. (77->43), mult. (161->55), div. (0->0), fcn. (137->6), ass. (0->24)
t11 = sin(qJ(1));
t14 = cos(qJ(1));
t34 = -g(1) * t11 + g(2) * t14;
t35 = -m(4) - m(5);
t33 = m(3) - t35;
t32 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t10 = sin(qJ(3));
t13 = cos(qJ(3));
t17 = mrSges(4,1) * t10 + mrSges(4,2) * t13;
t31 = mrSges(2,2) - mrSges(3,3) - m(5) * (pkin(3) * t10 - pkin(6) * t13) + t13 * mrSges(5,3) - t17;
t28 = t14 * pkin(1) + t11 * qJ(2);
t9 = sin(qJ(4));
t25 = t11 * t9;
t24 = t14 * t9;
t12 = cos(qJ(4));
t23 = t11 * t12;
t22 = t12 * t14;
t19 = m(5) * pkin(6) + mrSges(5,3);
t16 = m(5) * pkin(3) + mrSges(5,1) * t12 - mrSges(5,2) * t9;
t4 = t10 * t22 - t25;
t3 = t10 * t24 + t23;
t2 = t10 * t23 + t24;
t1 = -t10 * t25 + t22;
t5 = [(-m(3) * t28 - t2 * mrSges(5,1) - t1 * mrSges(5,2) + t35 * (t14 * pkin(5) + t28) + t32 * t14 + t31 * t11) * g(2) + (-t4 * mrSges(5,1) + t3 * mrSges(5,2) + (m(3) * pkin(1) + t35 * (-pkin(1) - pkin(5)) - t32) * t11 + (-t33 * qJ(2) + t31) * t14) * g(1), t34 * t33, (t16 * t10 - t19 * t13 + t17) * g(3) + ((mrSges(4,1) + t16) * t13 + (-mrSges(4,2) + t19) * t10) * t34, -g(1) * (mrSges(5,1) * t1 - mrSges(5,2) * t2) - g(2) * (mrSges(5,1) * t3 + mrSges(5,2) * t4) - g(3) * (-mrSges(5,1) * t9 - mrSges(5,2) * t12) * t13];
taug = t5(:);

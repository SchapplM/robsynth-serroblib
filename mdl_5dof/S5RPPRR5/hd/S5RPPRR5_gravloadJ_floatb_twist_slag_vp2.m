% Calculate Gravitation load on the joints for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:30
% DurationCPUTime: 0.25s
% Computational Cost: add. (201->47), mult. (189->56), div. (0->0), fcn. (184->8), ass. (0->26)
t12 = sin(qJ(5));
t14 = cos(qJ(5));
t20 = -mrSges(6,1) * t14 + mrSges(6,2) * t12;
t36 = m(6) * pkin(4) + mrSges(5,1) - t20;
t35 = mrSges(3,1) + mrSges(4,1);
t34 = mrSges(3,2) - mrSges(4,3);
t28 = m(4) + m(5) + m(6);
t27 = qJ(1) + pkin(8);
t10 = cos(t27);
t23 = sin(t27);
t29 = sin(qJ(4));
t30 = cos(qJ(4));
t1 = -t10 * t30 - t23 * t29;
t2 = t10 * t29 - t23 * t30;
t33 = t2 * mrSges(5,2) + t36 * t1;
t32 = t1 * mrSges(5,2) - t36 * t2;
t13 = sin(qJ(1));
t31 = t13 * pkin(1);
t15 = cos(qJ(1));
t11 = t15 * pkin(1);
t26 = t10 * pkin(2) + t23 * qJ(3) + t11;
t25 = -m(6) * pkin(7) - mrSges(6,3);
t24 = t10 * pkin(3) + t26;
t17 = -t23 * pkin(2) + t10 * qJ(3) - t31;
t16 = -t23 * pkin(3) + t17;
t3 = [(-mrSges(2,1) * t15 + t13 * mrSges(2,2) - m(3) * t11 - m(4) * t26 - m(5) * t24 - m(6) * (pkin(7) * t2 + t24) - t2 * mrSges(6,3) + t34 * t23 - t35 * t10 + t33) * g(2) + (t13 * mrSges(2,1) + mrSges(2,2) * t15 + m(3) * t31 - m(4) * t17 - m(5) * t16 - m(6) * (t1 * pkin(7) + t16) - t1 * mrSges(6,3) + t35 * t23 + t34 * t10 + t32) * g(1), (-m(3) - t28) * g(3), (-g(1) * t23 + g(2) * t10) * t28, (-t25 * t2 - t33) * g(2) + (-t25 * t1 - t32) * g(1), -g(3) * t20 + (-g(1) * t1 - g(2) * t2) * (mrSges(6,1) * t12 + mrSges(6,2) * t14)];
taug = t3(:);

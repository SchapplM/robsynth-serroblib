% Calculate Gravitation load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:19
% DurationCPUTime: 0.30s
% Computational Cost: add. (133->45), mult. (191->49), div. (0->0), fcn. (143->6), ass. (0->22)
t43 = mrSges(5,1) + mrSges(6,1);
t42 = mrSges(5,2) - mrSges(6,3);
t39 = m(5) + m(6);
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t9 = qJ(3) + pkin(7);
t4 = sin(t9);
t5 = cos(t9);
t41 = t11 * mrSges(4,1) + t13 * mrSges(4,2) + t4 * t43 + t42 * t5;
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t37 = -g(1) * t12 + g(2) * t14;
t40 = -m(3) - m(4);
t36 = mrSges(2,1) + mrSges(6,2) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t23 = pkin(4) * t4 - qJ(5) * t5;
t35 = -m(6) * t23 + mrSges(2,2) - mrSges(3,3) - t41;
t33 = t14 * pkin(1) + t12 * qJ(2);
t32 = pkin(3) * t11;
t7 = t14 * qJ(2);
t27 = -t12 * pkin(1) + t7;
t10 = -qJ(4) - pkin(6);
t1 = [(t40 * t33 - t39 * (-t10 * t14 + t12 * t32 + t33) + (-m(4) * pkin(6) - t36) * t14 + t35 * t12) * g(2) + (-m(3) * t27 - m(4) * t7 - t39 * (t12 * t10 + t14 * t32 + t27) + t35 * t14 + (-m(4) * (-pkin(1) - pkin(6)) + t36) * t12) * g(1), t37 * (t39 - t40), (m(5) * t32 - m(6) * (-t23 - t32) + t41) * g(3) + t37 * (-mrSges(4,2) * t11 + (m(6) * pkin(4) + t43) * t5 + (m(6) * qJ(5) - t42) * t4 + (pkin(3) * t39 + mrSges(4,1)) * t13), t39 * (-g(1) * t14 - g(2) * t12), (-g(3) * t4 - t37 * t5) * m(6)];
taug = t1(:);

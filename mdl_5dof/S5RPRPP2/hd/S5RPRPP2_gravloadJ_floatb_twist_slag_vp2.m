% Calculate Gravitation load on the joints for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:42
% DurationCPUTime: 0.34s
% Computational Cost: add. (172->50), mult. (192->51), div. (0->0), fcn. (147->6), ass. (0->23)
t48 = -mrSges(4,1) - mrSges(5,1);
t47 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t13 = qJ(1) + pkin(7);
t8 = sin(t13);
t9 = cos(t13);
t40 = g(1) * t9 + g(2) * t8;
t33 = m(5) + m(6);
t43 = -m(4) - t33;
t14 = sin(qJ(3));
t10 = t14 * qJ(4);
t16 = cos(qJ(3));
t11 = t16 * pkin(3);
t30 = t11 + t10;
t39 = t47 * t14 + t48 * t16;
t37 = -mrSges(3,1) + t39;
t36 = m(6) * qJ(5) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t15 = sin(qJ(1));
t32 = pkin(1) * t15;
t17 = cos(qJ(1));
t12 = t17 * pkin(1);
t26 = -pkin(2) - t10;
t23 = m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1);
t1 = [(-m(3) * t12 - mrSges(2,1) * t17 + t15 * mrSges(2,2) + t36 * t8 + t43 * (t9 * pkin(2) + t8 * pkin(6) + t12) + (-t33 * t30 - (m(6) * pkin(4) + mrSges(6,1)) * t16 + t37) * t9) * g(2) + (m(3) * t32 + t15 * mrSges(2,1) + mrSges(2,2) * t17 + t43 * (t9 * pkin(6) - t32) + t36 * t9 + (m(4) * pkin(2) - m(5) * (t26 - t11) - m(6) * t26 - t23 * t16 - t37) * t8) * g(1), (-m(3) + t43) * g(3), (-m(5) * t30 - m(6) * (pkin(4) * t16 + t30) - t16 * mrSges(6,1) + t39) * g(3) + ((m(5) * pkin(3) - t23 - t48) * t14 + (-qJ(4) * t33 + t47) * t16) * t40, (g(3) * t16 - t14 * t40) * t33, (g(1) * t8 - g(2) * t9) * m(6)];
taug = t1(:);

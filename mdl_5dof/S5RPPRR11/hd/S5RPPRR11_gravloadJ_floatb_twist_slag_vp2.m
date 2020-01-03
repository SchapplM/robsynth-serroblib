% Calculate Gravitation load on the joints for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:32
% EndTime: 2019-12-31 18:05:33
% DurationCPUTime: 0.27s
% Computational Cost: add. (98->52), mult. (195->63), div. (0->0), fcn. (159->6), ass. (0->26)
t38 = -m(5) - m(6);
t11 = sin(qJ(1));
t14 = cos(qJ(1));
t37 = -g(1) * t14 - g(2) * t11;
t25 = -m(4) + t38;
t36 = mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t10 = sin(qJ(4));
t13 = cos(qJ(4));
t17 = mrSges(5,1) * t10 + mrSges(5,2) * t13;
t35 = m(6) * (pkin(4) * t10 - pkin(7) * t13) + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t17 - t13 * mrSges(6,3);
t34 = t14 * pkin(1) + t11 * qJ(2);
t9 = sin(qJ(5));
t31 = t11 * t9;
t30 = t14 * t9;
t12 = cos(qJ(5));
t29 = t11 * t12;
t28 = t12 * t14;
t24 = t14 * qJ(3) + t34;
t22 = m(6) * pkin(7) + mrSges(6,3);
t16 = m(6) * pkin(4) + t12 * mrSges(6,1) - t9 * mrSges(6,2);
t7 = t14 * qJ(2);
t4 = t10 * t28 - t31;
t3 = -t10 * t30 - t29;
t2 = -t10 * t29 - t30;
t1 = t10 * t31 - t28;
t5 = [(-m(3) * t34 - m(4) * t24 - t4 * mrSges(6,1) - t3 * mrSges(6,2) + t38 * (-t11 * pkin(6) + t24) + t36 * t11 - t35 * t14) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (-m(3) - m(4)) * t7 + t38 * (-pkin(6) * t14 + t7) + t36 * t14 + (m(3) * pkin(1) + t25 * (-pkin(1) - qJ(3)) + t35) * t11) * g(1), (-g(1) * t11 + g(2) * t14) * (m(3) - t25), -t37 * t25, (t16 * t10 - t22 * t13 + t17) * g(3) + ((mrSges(5,1) + t16) * t13 + (-mrSges(5,2) + t22) * t10) * t37, -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t9 - mrSges(6,2) * t12) * t13];
taug = t5(:);

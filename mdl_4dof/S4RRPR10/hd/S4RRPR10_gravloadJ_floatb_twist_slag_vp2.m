% Calculate Gravitation load on the joints for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:11
% DurationCPUTime: 0.34s
% Computational Cost: add. (98->57), mult. (212->67), div. (0->0), fcn. (182->6), ass. (0->31)
t51 = -mrSges(3,1) + mrSges(4,2);
t50 = mrSges(3,2) - mrSges(4,3);
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t43 = g(1) * t19 + g(2) * t16;
t46 = m(4) + m(5);
t15 = sin(qJ(2));
t18 = cos(qJ(2));
t44 = t50 * t15 + t51 * t18;
t42 = -mrSges(2,1) + t44;
t40 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t37 = g(3) * t18;
t11 = t18 * pkin(2);
t9 = t15 * qJ(3);
t36 = t11 + t9;
t14 = sin(qJ(4));
t35 = t16 * t14;
t17 = cos(qJ(4));
t34 = t16 * t17;
t33 = t19 * t15;
t32 = t19 * t18;
t31 = t19 * pkin(1) + t16 * pkin(5);
t29 = pkin(2) * t32 + t19 * t9 + t31;
t28 = -pkin(1) - t9;
t27 = m(5) * (-pkin(2) - pkin(6)) - mrSges(5,3);
t23 = t14 * mrSges(5,1) + t17 * mrSges(5,2);
t4 = -t15 * t35 + t17 * t19;
t3 = t14 * t19 + t15 * t34;
t2 = t14 * t33 + t34;
t1 = t17 * t33 - t35;
t5 = [(-m(3) * t31 - m(4) * t29 - m(5) * (pkin(6) * t32 + t29) - t2 * mrSges(5,1) - t1 * mrSges(5,2) - mrSges(5,3) * t32 + t42 * t19 + t40 * t16) * g(2) + (-t4 * mrSges(5,1) + t3 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * (t28 - t11) - m(5) * t28 - t27 * t18 - t42) * t16 + ((-m(3) - t46) * pkin(5) + t40) * t19) * g(1), (-m(4) * t36 - m(5) * (pkin(6) * t18 + t36) - t18 * mrSges(5,3) - t23 * t15 + t44) * g(3) + ((m(4) * pkin(2) - t27 - t51) * t15 + (-qJ(3) * t46 - t23 + t50) * t18) * t43, (-t15 * t43 + t37) * t46, -g(1) * (mrSges(5,1) * t1 - mrSges(5,2) * t2) - g(2) * (mrSges(5,1) * t3 + mrSges(5,2) * t4) - (-mrSges(5,1) * t17 + mrSges(5,2) * t14) * t37];
taug = t5(:);

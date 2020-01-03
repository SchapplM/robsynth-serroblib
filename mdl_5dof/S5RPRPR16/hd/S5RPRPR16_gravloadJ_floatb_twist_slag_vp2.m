% Calculate Gravitation load on the joints for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR16_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:46
% DurationCPUTime: 0.43s
% Computational Cost: add. (119->72), mult. (246->85), div. (0->0), fcn. (204->6), ass. (0->33)
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t50 = t18 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * t15;
t49 = -m(3) - m(4);
t46 = m(5) + m(6);
t10 = t18 * qJ(4);
t48 = mrSges(2,2) - mrSges(3,3) - (-m(5) * qJ(4) - mrSges(5,3)) * t18 - m(6) * (t15 * pkin(7) - t10) - t15 * mrSges(6,3) - t50;
t47 = m(6) * pkin(4) + mrSges(2,1) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t19 = cos(qJ(1));
t40 = g(2) * t19;
t16 = sin(qJ(1));
t41 = g(1) * t16;
t45 = -t41 + t40;
t39 = g(3) * t15;
t38 = t15 * pkin(3);
t36 = t19 * pkin(1) + t16 * qJ(2);
t34 = t16 * t18;
t14 = sin(qJ(5));
t33 = t19 * t14;
t17 = cos(qJ(5));
t32 = t19 * t17;
t31 = qJ(4) * t15;
t30 = t19 * pkin(6) + t36;
t27 = m(6) * (-pkin(3) - pkin(7)) - mrSges(6,3);
t25 = mrSges(4,1) * t18 - mrSges(4,2) * t15;
t23 = t14 * mrSges(6,1) + t17 * mrSges(6,2);
t22 = -t18 * mrSges(5,2) + t15 * mrSges(5,3);
t11 = t19 * qJ(2);
t4 = -t14 * t34 + t32;
t3 = -t17 * t34 - t33;
t2 = -t16 * t17 - t18 * t33;
t1 = t16 * t14 - t18 * t32;
t5 = [(-m(3) * t36 - m(4) * t30 - t4 * mrSges(6,1) - t3 * mrSges(6,2) - t46 * (t16 * t38 + t30) - t47 * t19 + t48 * t16) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) - t46 * (t19 * t38 + t11) + t49 * t11 + (m(3) * pkin(1) + (-m(4) - t46) * (-pkin(1) - pkin(6)) + t47) * t16 + t48 * t19) * g(1), t45 * (t46 - t49), -t25 * t41 + (-t46 * (pkin(3) * t34 + t16 * t31) + (-t22 - (m(6) * pkin(7) + mrSges(6,3)) * t18 - t23 * t15) * t16) * g(1) + (t25 - m(5) * (-pkin(3) * t18 - t31) + t22 - t27 * t18 - (-m(6) * qJ(4) - t23) * t15) * t40 + (-m(5) * (t10 - t38) - m(6) * t10 - t27 * t15 + (-mrSges(5,3) - t23) * t18 + t50) * g(3), (-t18 * t45 - t39) * t46, -g(1) * (t3 * mrSges(6,1) - t4 * mrSges(6,2)) - g(2) * (-t1 * mrSges(6,1) + t2 * mrSges(6,2)) - (mrSges(6,1) * t17 - mrSges(6,2) * t14) * t39];
taug = t5(:);

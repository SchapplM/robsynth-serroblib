% Calculate Gravitation load on the joints for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:42
% DurationCPUTime: 0.22s
% Computational Cost: add. (201->46), mult. (131->48), div. (0->0), fcn. (92->8), ass. (0->28)
t51 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t49 = mrSges(6,1) * t22 + mrSges(6,2) * t24;
t50 = -t49 + mrSges(4,2) - mrSges(5,3);
t48 = m(5) + m(6);
t21 = qJ(1) + pkin(8);
t19 = qJ(3) + t21;
t15 = sin(t19);
t16 = cos(t19);
t47 = -g(1) * t15 + g(2) * t16;
t46 = t50 * t16 + (-m(6) * (-pkin(3) - pkin(7)) - t51) * t15;
t45 = t50 * t15 + t51 * t16;
t23 = sin(qJ(1));
t44 = pkin(1) * t23;
t43 = pkin(3) * t15;
t25 = cos(qJ(1));
t20 = t25 * pkin(1);
t40 = t16 * pkin(3) + t15 * qJ(4);
t18 = cos(t21);
t36 = pkin(2) * t18 + t20;
t32 = t36 + t40;
t17 = sin(t21);
t31 = -pkin(2) * t17 - t44;
t6 = t16 * qJ(4);
t27 = t31 + t6;
t12 = t16 * pkin(7);
t1 = [(-mrSges(2,1) * t25 + t23 * mrSges(2,2) - m(3) * t20 - t18 * mrSges(3,1) + t17 * mrSges(3,2) - m(4) * t36 - m(5) * t32 - m(6) * (t12 + t32) + t45) * g(2) + (t23 * mrSges(2,1) + mrSges(2,2) * t25 + m(3) * t44 + mrSges(3,1) * t17 + mrSges(3,2) * t18 - m(4) * t31 - m(5) * (t27 - t43) - m(6) * t27 + t46) * g(1), (-m(3) - m(4) - t48) * g(3), (-m(5) * t40 - m(6) * (t12 + t40) + t45) * g(2) + (-m(5) * (t6 - t43) - m(6) * t6 + t46) * g(1), t48 * t47, g(3) * t49 + t47 * (mrSges(6,1) * t24 - mrSges(6,2) * t22)];
taug = t1(:);

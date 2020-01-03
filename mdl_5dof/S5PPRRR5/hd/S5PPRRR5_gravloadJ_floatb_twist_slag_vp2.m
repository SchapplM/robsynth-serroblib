% Calculate Gravitation load on the joints for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:38
% EndTime: 2019-12-31 17:35:39
% DurationCPUTime: 0.18s
% Computational Cost: add. (154->32), mult. (178->44), div. (0->0), fcn. (182->8), ass. (0->24)
t21 = sin(qJ(5));
t23 = cos(qJ(5));
t44 = -mrSges(6,1) * t23 + t21 * mrSges(6,2);
t46 = mrSges(5,1) - t44;
t45 = -mrSges(5,2) + mrSges(6,3);
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t14 = -t19 * t24 + t20 * t22;
t35 = qJ(3) + qJ(4);
t30 = sin(t35);
t31 = cos(t35);
t11 = -t19 * t30 - t20 * t31;
t12 = -t19 * t31 + t20 * t30;
t43 = t45 * t11 + t46 * t12;
t42 = -t46 * t11 + t45 * t12;
t34 = m(3) + m(4) + m(5) + m(6);
t33 = t11 * pkin(4) - t12 * pkin(7);
t32 = -t12 * pkin(4) - pkin(7) * t11;
t29 = t14 * pkin(3);
t13 = -t19 * t22 - t20 * t24;
t25 = t13 * pkin(3);
t1 = [(-m(2) - t34) * g(3), (-g(1) * t19 + g(2) * t20) * t34, (-mrSges(4,1) * t13 - mrSges(4,2) * t14 - m(5) * t25 - m(6) * (t25 + t33) + t42) * g(2) + (mrSges(4,1) * t14 - mrSges(4,2) * t13 + m(5) * t29 - m(6) * (-t29 + t32) + t43) * g(1), (-m(6) * t33 + t42) * g(2) + (-m(6) * t32 + t43) * g(1), -g(3) * t44 + (-g(1) * t11 - g(2) * t12) * (mrSges(6,1) * t21 + mrSges(6,2) * t23)];
taug = t1(:);

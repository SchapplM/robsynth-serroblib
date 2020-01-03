% Calculate Gravitation load on the joints for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:04
% DurationCPUTime: 0.38s
% Computational Cost: add. (107->41), mult. (254->56), div. (0->0), fcn. (243->8), ass. (0->21)
t45 = m(5) + m(6);
t54 = -mrSges(3,1) - mrSges(4,1);
t53 = mrSges(3,2) - mrSges(4,3);
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t42 = g(1) * t17 + g(2) * t16;
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t32 = sin(pkin(8));
t33 = cos(pkin(8));
t22 = -t19 * t33 + t21 * t32;
t50 = t22 * g(3);
t31 = m(4) + t45;
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t44 = -m(6) * pkin(4) - t20 * mrSges(6,1) + t18 * mrSges(6,2) - mrSges(5,1);
t43 = t19 * t32 + t21 * t33;
t35 = t21 * pkin(2) + t19 * qJ(3);
t3 = t43 * t17;
t1 = t43 * t16;
t2 = [(-m(2) - m(3) - t31) * g(3), (-m(4) * t35 - t45 * (t21 * pkin(3) + t35) + t44 * t43 + t54 * t21 + t53 * t19) * g(3) + (g(1) * t3 + g(2) * t1 - t50) * (m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3)) + (t22 * t44 + (-qJ(3) * t31 + t53) * t21 + (-t45 * (-pkin(2) - pkin(3)) + m(4) * pkin(2) - t54) * t19) * t42, (g(3) * t21 - t42 * t19) * t31, t45 * (g(1) * t16 - g(2) * t17), -g(1) * ((-t16 * t20 - t18 * t3) * mrSges(6,1) + (t16 * t18 - t20 * t3) * mrSges(6,2)) - g(2) * ((-t1 * t18 + t17 * t20) * mrSges(6,1) + (-t1 * t20 - t17 * t18) * mrSges(6,2)) + (-t18 * mrSges(6,1) - t20 * mrSges(6,2)) * t50];
taug = t2(:);

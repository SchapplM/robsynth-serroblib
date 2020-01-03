% Calculate Gravitation load on the joints for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:41
% EndTime: 2019-12-31 18:40:42
% DurationCPUTime: 0.28s
% Computational Cost: add. (248->53), mult. (181->55), div. (0->0), fcn. (137->8), ass. (0->30)
t61 = mrSges(5,1) + mrSges(6,1);
t60 = -mrSges(5,2) + mrSges(6,3);
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t59 = t60 * t23 + t61 * t25;
t32 = pkin(4) * t25 + qJ(5) * t23;
t58 = -mrSges(4,1) - t59;
t57 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t22 = qJ(1) + pkin(8);
t20 = qJ(3) + t22;
t17 = cos(t20);
t55 = t32 * t17;
t16 = sin(t20);
t54 = g(1) * t17 + g(2) * t16;
t53 = t57 * t17 + (-m(6) * (-pkin(3) - t32) - t58) * t16;
t52 = t57 * t16 + t58 * t17;
t24 = sin(qJ(1));
t51 = pkin(1) * t24;
t50 = pkin(3) * t16;
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t42 = t17 * pkin(3) + t16 * pkin(7);
t19 = cos(t22);
t41 = pkin(2) * t19 + t21;
t37 = t41 + t42;
t18 = sin(t22);
t36 = -pkin(2) * t18 - t51;
t13 = t17 * pkin(7);
t29 = t13 + t36;
t1 = [(-mrSges(2,1) * t26 + t24 * mrSges(2,2) - m(3) * t21 - t19 * mrSges(3,1) + t18 * mrSges(3,2) - m(4) * t41 - m(5) * t37 - m(6) * (t37 + t55) + t52) * g(2) + (t24 * mrSges(2,1) + mrSges(2,2) * t26 + m(3) * t51 + mrSges(3,1) * t18 + mrSges(3,2) * t19 - m(4) * t36 - m(5) * (t29 - t50) - m(6) * t29 + t53) * g(1), (-m(3) - m(4) - m(5) - m(6)) * g(3), (-m(5) * t42 - m(6) * (t42 + t55) + t52) * g(2) + (-m(5) * (t13 - t50) - m(6) * t13 + t53) * g(1), (-m(6) * t32 - t59) * g(3) + ((-m(6) * qJ(5) - t60) * t25 + (m(6) * pkin(4) + t61) * t23) * t54, (g(3) * t25 - t23 * t54) * m(6)];
taug = t1(:);

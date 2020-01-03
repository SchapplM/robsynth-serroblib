% Calculate Gravitation load on the joints for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:09
% DurationCPUTime: 0.33s
% Computational Cost: add. (180->43), mult. (341->50), div. (0->0), fcn. (363->6), ass. (0->26)
t14 = sin(qJ(4));
t15 = cos(qJ(4));
t59 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t60 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t53 = t59 * t14 + t60 * t15;
t63 = mrSges(4,1) + t53;
t33 = sin(qJ(3));
t34 = sin(qJ(1));
t35 = cos(qJ(3));
t36 = cos(qJ(1));
t5 = -t34 * t33 - t36 * t35;
t58 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3);
t6 = t36 * t33 - t34 * t35;
t62 = -t58 * t5 - t63 * t6;
t61 = t63 * t5 - t58 * t6;
t57 = -m(5) - m(6);
t56 = mrSges(2,1) + mrSges(3,1);
t55 = mrSges(2,2) - mrSges(3,3);
t52 = -t6 * pkin(3) - t5 * pkin(7);
t51 = t5 * pkin(3) - t6 * pkin(7);
t50 = g(1) * t5 + g(2) * t6;
t32 = t36 * pkin(1) + t34 * qJ(2);
t31 = t36 * pkin(2) + t32;
t28 = -t34 * pkin(1) + t36 * qJ(2);
t20 = -t34 * pkin(2) + t28;
t1 = [(-m(3) * t32 - m(4) * t31 - t56 * t36 + t55 * t34 + t57 * (t31 - t51) + t61) * g(2) + (-m(3) * t28 - m(4) * t20 + t55 * t36 + t56 * t34 + t57 * (t20 - t52) + t62) * g(1), (-t34 * g(1) + t36 * g(2)) * (m(3) + m(4) - t57), (t57 * t51 - t61) * g(2) + (t57 * t52 - t62) * g(1), t53 * g(3) + (-t60 * t14 + t59 * t15) * t50, (-g(3) * t15 + t14 * t50) * m(6)];
taug = t1(:);

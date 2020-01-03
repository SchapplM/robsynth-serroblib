% Calculate Gravitation load on the joints for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:40
% EndTime: 2019-12-31 21:59:42
% DurationCPUTime: 0.75s
% Computational Cost: add. (308->94), mult. (454->113), div. (0->0), fcn. (420->8), ass. (0->51)
t86 = mrSges(5,1) + mrSges(6,1);
t85 = mrSges(5,2) + mrSges(6,2);
t87 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t29 = qJ(3) + qJ(4);
t23 = cos(t29);
t33 = cos(qJ(3));
t25 = t33 * pkin(3);
t19 = pkin(4) * t23 + t25;
t17 = pkin(2) + t19;
t21 = t25 + pkin(2);
t22 = sin(t29);
t30 = sin(qJ(3));
t84 = -m(4) * pkin(2) - m(5) * t21 - m(6) * t17 - t33 * mrSges(4,1) + t30 * mrSges(4,2) + t85 * t22 - t86 * t23;
t36 = -pkin(8) - pkin(7);
t28 = -qJ(5) + t36;
t83 = -m(4) * pkin(7) + m(5) * t36 + m(6) * t28 - t87;
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t82 = g(1) * t35 + g(2) * t32;
t70 = m(5) * pkin(3);
t64 = pkin(3) * t30;
t18 = pkin(4) * t22 + t64;
t81 = m(6) * t18;
t80 = mrSges(4,1) + t70;
t77 = mrSges(5,1) * t22 + t85 * t23;
t34 = cos(qJ(2));
t53 = t34 * t35;
t11 = -t22 * t53 + t23 * t32;
t12 = t22 * t32 + t23 * t53;
t76 = -t86 * t11 + t85 * t12;
t54 = t32 * t34;
t10 = t22 * t35 - t23 * t54;
t9 = t22 * t54 + t23 * t35;
t75 = -t85 * t10 + t86 * t9;
t74 = -m(3) - m(4) - m(5) - m(6);
t72 = mrSges(2,2) - mrSges(3,3) - t81;
t31 = sin(qJ(2));
t47 = t34 * mrSges(3,1) - t31 * mrSges(3,2);
t71 = t87 * t31 + mrSges(2,1) + t47;
t69 = m(6) * pkin(4);
t61 = g(3) * t31;
t59 = t30 * t32;
t58 = t30 * t35;
t48 = pkin(2) * t34 + pkin(7) * t31;
t44 = t17 * t34 - t28 * t31;
t43 = t34 * t21 - t31 * t36;
t15 = -t30 * t53 + t32 * t33;
t13 = t30 * t54 + t33 * t35;
t16 = t33 * t53 + t59;
t14 = -t33 * t54 + t58;
t1 = [(-t59 * t70 - t16 * mrSges(4,1) - t15 * mrSges(4,2) + t74 * (t35 * pkin(1) + t32 * pkin(6)) + t72 * t32 - t86 * t12 - t85 * t11 + (-m(4) * t48 - m(5) * t43 - m(6) * t44 - t71) * t35) * g(2) + (-t58 * t70 - t14 * mrSges(4,1) - t13 * mrSges(4,2) - t85 * t9 - t86 * t10 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t48) - m(5) * (-pkin(1) - t43) - m(6) * (-pkin(1) - t44) + t71) * t32 + (t74 * pkin(6) + t72) * t35) * g(1), -g(3) * t47 + (t84 * g(3) + t82 * (mrSges(3,2) + t83)) * t34 + (t83 * g(3) + t82 * (mrSges(3,1) - t84)) * t31, (m(5) * t64 + mrSges(4,1) * t30 + mrSges(6,1) * t22 + mrSges(4,2) * t33 + t77 + t81) * t61 + (-mrSges(4,2) * t14 - m(6) * (-t18 * t54 - t19 * t35) + t80 * t13 + t75) * g(2) + (mrSges(4,2) * t16 - m(6) * (-t18 * t53 + t19 * t32) - t80 * t15 + t76) * g(1), (-(-mrSges(6,1) - t69) * t22 + t77) * t61 + (t9 * t69 + t75) * g(2) + (-t11 * t69 + t76) * g(1), (g(3) * t34 - t82 * t31) * m(6)];
taug = t1(:);

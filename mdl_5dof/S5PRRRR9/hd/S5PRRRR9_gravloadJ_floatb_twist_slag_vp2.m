% Calculate Gravitation load on the joints for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:03
% EndTime: 2019-12-05 17:19:07
% DurationCPUTime: 0.80s
% Computational Cost: add. (371->73), mult. (818->110), div. (0->0), fcn. (956->12), ass. (0->38)
t27 = qJ(4) + qJ(5);
t25 = sin(t27);
t26 = cos(t27);
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t69 = -mrSges(4,1) - m(6) * (pkin(4) * t33 + pkin(3)) - mrSges(6,1) * t26 + mrSges(6,2) * t25 - m(5) * pkin(3) - t33 * mrSges(5,1) + t30 * mrSges(5,2);
t68 = mrSges(4,2) + m(6) * (-pkin(9) - pkin(8)) - mrSges(6,3) - m(5) * pkin(8) - mrSges(5,3);
t74 = -m(6) * pkin(4) - mrSges(5,1);
t31 = sin(qJ(3));
t34 = cos(qJ(3));
t83 = t68 * t31 + t69 * t34 - mrSges(3,1);
t70 = m(4) + m(5) + m(6);
t78 = -t25 * mrSges(6,1) - t33 * mrSges(5,2) - t26 * mrSges(6,2) + t74 * t30 + mrSges(3,2) - mrSges(4,3);
t77 = pkin(2) * t70 - t83;
t65 = -t70 * pkin(7) + t78;
t28 = sin(pkin(10));
t32 = sin(qJ(2));
t35 = cos(qJ(2));
t54 = cos(pkin(10));
t55 = cos(pkin(5));
t44 = t55 * t54;
t13 = t28 * t32 - t35 * t44;
t14 = t28 * t35 + t32 * t44;
t29 = sin(pkin(5));
t50 = t29 * t54;
t8 = t14 * t34 - t31 * t50;
t63 = (t13 * t26 - t25 * t8) * mrSges(6,1) + (-t13 * t25 - t26 * t8) * mrSges(6,2);
t51 = t28 * t55;
t16 = -t32 * t51 + t54 * t35;
t10 = t28 * t29 * t31 + t16 * t34;
t15 = t54 * t32 + t35 * t51;
t62 = (-t10 * t25 + t15 * t26) * mrSges(6,1) + (-t10 * t26 - t15 * t25) * mrSges(6,2);
t58 = t29 * t34;
t18 = t55 * t31 + t32 * t58;
t57 = t29 * t35;
t61 = (-t18 * t25 - t26 * t57) * mrSges(6,1) + (-t18 * t26 + t25 * t57) * mrSges(6,2);
t59 = t29 * t32;
t1 = [(-m(2) - m(3) - t70) * g(3), (-t70 * (pkin(2) * t57 + pkin(7) * t59) + (t78 * t32 + t83 * t35) * t29) * g(3) + (t77 * t13 + t65 * t14) * g(2) + (t77 * t15 + t65 * t16) * g(1), (t68 * t18 + t69 * (-t31 * t59 + t55 * t34)) * g(3) + (t68 * t8 + t69 * (-t14 * t31 - t34 * t50)) * g(2) + (t69 * (-t16 * t31 + t28 * t58) + t68 * t10) * g(1), (-(-t18 * t33 + t30 * t57) * mrSges(5,2) - t61 + t74 * (-t18 * t30 - t33 * t57)) * g(3) + (-(-t13 * t30 - t33 * t8) * mrSges(5,2) - t63 + t74 * (t13 * t33 - t30 * t8)) * g(2) + (-(-t10 * t33 - t15 * t30) * mrSges(5,2) - t62 + t74 * (-t10 * t30 + t15 * t33)) * g(1), -g(1) * t62 - g(2) * t63 - g(3) * t61];
taug = t1(:);

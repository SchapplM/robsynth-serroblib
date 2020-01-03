% Calculate Gravitation load on the joints for
% S5RPRRP4
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:15
% EndTime: 2020-01-03 11:49:18
% DurationCPUTime: 0.55s
% Computational Cost: add. (232->77), mult. (337->95), div. (0->0), fcn. (316->8), ass. (0->45)
t73 = -mrSges(5,1) - mrSges(6,1);
t71 = mrSges(5,2) + mrSges(6,2);
t65 = m(5) * pkin(3);
t74 = -mrSges(4,1) - t65;
t72 = mrSges(2,2) - mrSges(3,3);
t31 = qJ(3) + qJ(4);
t24 = sin(t31);
t25 = cos(t31);
t70 = mrSges(5,1) * t24 + t71 * t25;
t33 = cos(pkin(8));
t37 = cos(qJ(1));
t50 = t37 * t24;
t35 = sin(qJ(1));
t54 = t35 * t25;
t11 = t33 * t50 - t54;
t49 = t37 * t25;
t55 = t35 * t24;
t12 = t33 * t49 + t55;
t69 = t73 * t11 - t71 * t12;
t10 = t33 * t54 - t50;
t9 = -t33 * t55 - t49;
t68 = t71 * t10 + t73 * t9;
t67 = m(3) + m(4) + m(5) + m(6);
t36 = cos(qJ(3));
t28 = t36 * pkin(3);
t19 = pkin(4) * t25 + t28;
t32 = sin(pkin(8));
t38 = -pkin(7) - pkin(6);
t66 = -mrSges(2,1) + (-mrSges(3,1) - m(4) * pkin(2) - m(5) * (t28 + pkin(2)) - m(6) * (pkin(2) + t19)) * t33 + (mrSges(3,2) - m(4) * pkin(6) - mrSges(4,3) + m(5) * t38 - mrSges(5,3) - m(6) * (qJ(5) - t38) - mrSges(6,3)) * t32;
t64 = m(6) * pkin(4);
t59 = g(1) * t32;
t34 = sin(qJ(3));
t58 = t34 * pkin(3);
t18 = pkin(4) * t24 + t58;
t56 = t35 * t18;
t53 = t35 * t34;
t52 = t35 * t36;
t51 = t37 * t18;
t48 = t37 * t34;
t47 = t37 * t36;
t15 = t33 * t48 - t52;
t13 = -t33 * t53 - t47;
t16 = t33 * t47 + t53;
t14 = t33 * t52 - t48;
t1 = [(t48 * t65 + m(6) * t51 - t14 * mrSges(4,1) - t13 * mrSges(4,2) - t71 * t9 - t67 * (t35 * pkin(1) - t37 * qJ(2)) - t72 * t37 + t73 * t10 + t66 * t35) * g(3) + (-t53 * t65 - m(6) * t56 - t16 * mrSges(4,1) + t15 * mrSges(4,2) - t67 * (t37 * pkin(1) + t35 * qJ(2)) + t72 * t35 + t73 * t12 + t71 * t11 + t66 * t37) * g(2), (g(2) * t37 + g(3) * t35) * t67, (m(5) * t58 + m(6) * t18 + mrSges(4,1) * t34 + mrSges(6,1) * t24 + mrSges(4,2) * t36 + t70) * t59 + (-t16 * mrSges(4,2) - m(6) * (-t35 * t19 + t33 * t51) + t74 * t15 + t69) * g(3) + (t14 * mrSges(4,2) - m(6) * (-t37 * t19 - t33 * t56) + t74 * t13 + t68) * g(2), (-(-mrSges(6,1) - t64) * t24 + t70) * t59 + (-t11 * t64 + t69) * g(3) + (-t9 * t64 + t68) * g(2), (g(1) * t33 + (-g(2) * t35 + g(3) * t37) * t32) * m(6)];
taug = t1(:);

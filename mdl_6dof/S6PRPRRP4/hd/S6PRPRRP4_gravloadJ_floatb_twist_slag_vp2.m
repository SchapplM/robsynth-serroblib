% Calculate Gravitation load on the joints for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:16
% EndTime: 2019-03-08 20:09:18
% DurationCPUTime: 0.72s
% Computational Cost: add. (579->86), mult. (1014->124), div. (0->0), fcn. (1176->12), ass. (0->53)
t55 = sin(qJ(5));
t65 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t109 = t65 * t55;
t49 = pkin(11) + qJ(4);
t48 = cos(t49);
t66 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t108 = t48 * t66;
t107 = mrSges(6,3) + mrSges(7,2);
t106 = -m(4) * qJ(3) - t66 * t55 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t100 = -m(6) - m(7);
t104 = -m(5) + t100;
t47 = sin(t49);
t53 = cos(pkin(11));
t103 = -t48 * mrSges(5,1) + t47 * mrSges(5,2) - mrSges(3,1) - m(4) * pkin(2) - t53 * mrSges(4,1) + sin(pkin(11)) * mrSges(4,2);
t57 = cos(qJ(5));
t102 = -t65 * t57 + t106;
t92 = pkin(4) * t48;
t101 = t100 * (-pkin(9) * t47 - t92) + t57 * t108 - t48 * t109 - t103 + t107 * t47;
t96 = -t66 * t57 - mrSges(5,1) + t109;
t94 = mrSges(5,2) - t107;
t51 = sin(pkin(10));
t52 = sin(pkin(6));
t85 = t51 * t52;
t56 = sin(qJ(2));
t84 = t52 * t56;
t58 = cos(qJ(2));
t83 = t52 * t58;
t82 = t57 * t58;
t79 = cos(pkin(6));
t78 = cos(pkin(10));
t77 = t47 * t83;
t76 = t55 * t83;
t75 = m(4) - t104;
t70 = t51 * t79;
t69 = t52 * t78;
t62 = t79 * t78;
t54 = -pkin(8) - qJ(3);
t46 = pkin(3) * t53 + pkin(2);
t38 = t46 * t83;
t37 = -t56 * t70 + t58 * t78;
t36 = t56 * t78 + t58 * t70;
t35 = t51 * t58 + t56 * t62;
t34 = t51 * t56 - t58 * t62;
t27 = t47 * t79 + t48 * t84;
t26 = -t47 * t84 + t48 * t79;
t15 = t27 * t55 + t52 * t82;
t14 = t37 * t48 + t47 * t85;
t13 = -t37 * t47 + t48 * t85;
t12 = t35 * t48 - t47 * t69;
t11 = -t35 * t47 - t48 * t69;
t3 = t14 * t55 - t36 * t57;
t1 = t12 * t55 - t34 * t57;
t2 = [(-m(2) - m(3) - t75) * g(3) (t104 * (-t34 * t46 - t35 * t54) + t102 * t35 + t101 * t34) * g(2) + (t104 * (-t36 * t46 - t37 * t54) + t102 * t37 + t101 * t36) * g(1) + (-m(5) * t38 - t107 * t77 + t100 * (pkin(9) * t77 - t54 * t84 + t83 * t92 + t38) + t65 * (t48 * t76 - t57 * t84) + (-t82 * t108 + t103 * t58 + (m(5) * t54 + t106) * t56) * t52) * g(3) (-g(1) * t36 - g(2) * t34 + g(3) * t83) * t75 (t100 * (t26 * pkin(4) + pkin(9) * t27) + t94 * t27 + t96 * t26) * g(3) + (t100 * (t11 * pkin(4) + pkin(9) * t12) + t94 * t12 + t96 * t11) * g(2) + (t100 * (t13 * pkin(4) + pkin(9) * t14) + t94 * t14 + t96 * t13) * g(1) (t65 * (t27 * t57 - t76) + t66 * t15) * g(3) + (t65 * (t12 * t57 + t34 * t55) + t66 * t1) * g(2) + (t65 * (t14 * t57 + t36 * t55) + t66 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);

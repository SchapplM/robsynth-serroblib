% Calculate Gravitation load on the joints for
% S6PRPRRP3
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:00:41
% EndTime: 2018-11-23 15:00:42
% DurationCPUTime: 0.77s
% Computational Cost: add. (1040->90), mult. (1072->123), div. (0->0), fcn. (1020->16), ass. (0->53)
t54 = cos(qJ(5));
t104 = m(6) * pkin(4) + m(7) * (pkin(5) * t54 + pkin(4)) + mrSges(5,1);
t103 = m(6) * pkin(9) - m(7) * (-qJ(6) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t102 = mrSges(6,1) + mrSges(7,1);
t100 = -mrSges(6,2) - mrSges(7,2);
t46 = pkin(11) + qJ(4);
t44 = sin(t46);
t45 = cos(t46);
t49 = cos(pkin(11));
t93 = m(4) * pkin(2) + t49 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2) + mrSges(3,1) + t104 * t45 + t103 * t44;
t92 = m(7) * pkin(5);
t99 = -m(5) - m(6) - m(7);
t52 = sin(qJ(5));
t97 = t100 * t52 + t102 * t54 + t104;
t96 = -t92 - t102;
t95 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t77 = pkin(6) + qJ(2);
t67 = sin(t77) / 0.2e1;
t78 = pkin(6) - qJ(2);
t70 = sin(t78);
t34 = t67 - t70 / 0.2e1;
t48 = sin(pkin(10));
t55 = cos(qJ(2));
t80 = cos(pkin(10));
t63 = t80 * t34 + t48 * t55;
t91 = t63 * t52;
t62 = -t48 * t34 + t80 * t55;
t90 = t62 * t52;
t68 = cos(t77) / 0.2e1;
t71 = cos(t78);
t35 = t68 - t71 / 0.2e1;
t89 = t35 * t52;
t86 = t45 * t52;
t85 = t45 * t54;
t81 = cos(pkin(6));
t79 = sin(pkin(6));
t76 = m(4) - t99;
t72 = t48 * t79;
t64 = t80 * t79;
t56 = t71 / 0.2e1 + t68;
t53 = sin(qJ(2));
t51 = -pkin(8) - qJ(3);
t42 = pkin(3) * t49 + pkin(2);
t33 = t67 + t70 / 0.2e1;
t26 = t48 * t56 + t80 * t53;
t23 = t48 * t53 - t80 * t56;
t20 = -t35 * t45 + t81 * t44;
t19 = -t35 * t44 - t81 * t45;
t14 = t44 * t72 + t45 * t62;
t13 = t44 * t62 - t45 * t72;
t12 = -t44 * t64 + t45 * t63;
t11 = t44 * t63 + t45 * t64;
t1 = [(-m(2) - m(3) - t76) * g(3) (t89 * t92 + t99 * (t33 * t42 + t35 * t51) - t102 * (t33 * t85 - t89) + t100 * (-t33 * t86 - t35 * t54) - t95 * t35 - t93 * t33) * g(3) + (-t91 * t92 + t99 * (-t23 * t42 - t51 * t63) - t102 * (-t23 * t85 + t91) + t100 * (t23 * t86 + t54 * t63) + t95 * t63 + t93 * t23) * g(2) + (-t90 * t92 + t99 * (-t26 * t42 - t51 * t62) - t102 * (-t26 * t85 + t90) + t100 * (t26 * t86 + t54 * t62) + t95 * t62 + t93 * t26) * g(1) (-g(1) * t26 - g(2) * t23 + g(3) * t33) * t76 (-t103 * t20 + t97 * t19) * g(3) + (-t103 * t12 + t97 * t11) * g(2) + (-t103 * t14 + t97 * t13) * g(1) (t100 * (-t20 * t54 + t33 * t52) + t96 * (-t20 * t52 - t33 * t54)) * g(3) + (t100 * (-t12 * t54 - t23 * t52) + t96 * (-t12 * t52 + t23 * t54)) * g(2) + (t100 * (-t14 * t54 - t26 * t52) + t96 * (-t14 * t52 + t26 * t54)) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t19) * m(7)];
taug  = t1(:);

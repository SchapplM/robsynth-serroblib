% Calculate Gravitation load on the joints for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2018-11-23 15:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:12:19
% EndTime: 2018-11-23 15:12:20
% DurationCPUTime: 0.85s
% Computational Cost: add. (1271->92), mult. (1439->120), div. (0->0), fcn. (1420->16), ass. (0->60)
t122 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t59 = pkin(11) + qJ(5);
t57 = sin(t59);
t58 = cos(t59);
t60 = sin(pkin(11));
t61 = cos(pkin(11));
t84 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t113 = m(5) * pkin(3) + t61 * mrSges(5,1) - t60 * mrSges(5,2) + t122 * t57 + t84 * t58 + mrSges(4,1);
t112 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t116 = -m(6) - m(7);
t56 = pkin(4) * t61 + pkin(3);
t62 = -pkin(9) - qJ(4);
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t117 = t116 * (-t56 * t65 + t62 * t63) - t112 * t63 + mrSges(3,1) + t113 * t65;
t121 = -m(4) + t116;
t118 = -m(5) * pkin(8) - t61 * mrSges(5,2) - t84 * t57 + t122 * t58 + mrSges(3,2) - mrSges(4,3) + (pkin(4) * t116 - mrSges(5,1)) * t60;
t115 = m(5) - t116;
t94 = pkin(6) + qJ(2);
t85 = sin(t94);
t52 = t85 / 0.2e1;
t95 = pkin(6) - qJ(2);
t86 = sin(t95);
t100 = t52 - t86 / 0.2e1;
t99 = cos(pkin(6));
t98 = cos(pkin(10));
t97 = sin(pkin(6));
t96 = sin(pkin(10));
t66 = cos(qJ(2));
t89 = t98 * t66;
t88 = t96 * t66;
t87 = cos(t95);
t82 = cos(t94) / 0.2e1;
t81 = t86 / 0.2e1;
t78 = t98 * t97;
t77 = t97 * t96;
t70 = t87 / 0.2e1 + t82;
t69 = t81 - t85 / 0.2e1;
t64 = sin(qJ(2));
t49 = t82 - t87 / 0.2e1;
t48 = t52 + t81;
t47 = t48 * pkin(2);
t40 = -t49 * t65 + t99 * t63;
t39 = -t49 * t63 - t99 * t65;
t38 = t96 * t69 + t89;
t37 = -t96 * t100 + t89;
t36 = t98 * t64 + t96 * t70;
t35 = -t98 * t69 + t88;
t34 = t98 * t100 + t88;
t33 = t96 * t64 - t98 * t70;
t32 = t36 * pkin(2);
t31 = t33 * pkin(2);
t20 = t37 * t65 + t63 * t77;
t19 = t37 * t63 - t65 * t77;
t18 = t34 * t65 - t63 * t78;
t17 = t34 * t63 + t65 * t78;
t9 = t40 * t57 + t48 * t58;
t3 = t20 * t57 - t36 * t58;
t1 = t18 * t57 - t33 * t58;
t2 = [(-m(2) - m(3) - m(4) - t115) * g(3) (-m(5) * t47 + t121 * (-pkin(8) * t49 + t47) - t118 * t49 - t117 * t48) * g(3) + (m(5) * t31 + t121 * (pkin(8) * t35 - t31) + t118 * t35 + t117 * t33) * g(2) + (m(5) * t32 + t121 * (pkin(8) * t38 - t32) + t118 * t38 + t117 * t36) * g(1) (t116 * (-t39 * t56 - t40 * t62) + t112 * t40 + t113 * t39) * g(3) + (t116 * (-t17 * t56 - t18 * t62) + t112 * t18 + t113 * t17) * g(2) + (t116 * (-t19 * t56 - t20 * t62) + t112 * t20 + t113 * t19) * g(1), t115 * (-g(1) * t19 - g(2) * t17 - g(3) * t39) (t84 * t9 - t122 * (t40 * t58 - t48 * t57)) * g(3) + (-t122 * (t18 * t58 + t33 * t57) + t84 * t1) * g(2) + (-t122 * (t20 * t58 + t36 * t57) + t84 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);

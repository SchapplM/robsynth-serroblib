% Calculate Gravitation load on the joints for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:39
% EndTime: 2019-12-31 21:26:42
% DurationCPUTime: 0.93s
% Computational Cost: add. (475->110), mult. (878->158), div. (0->0), fcn. (1000->12), ass. (0->50)
t45 = sin(qJ(5));
t49 = cos(qJ(5));
t95 = m(6) * pkin(4) + t49 * mrSges(6,1) - t45 * mrSges(6,2) + mrSges(5,1);
t93 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t59 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t92 = t45 * mrSges(6,1) + t49 * mrSges(6,2) - t59;
t97 = m(5) + m(6);
t46 = sin(qJ(3));
t50 = cos(qJ(3));
t77 = cos(pkin(5));
t43 = sin(pkin(5));
t47 = sin(qJ(2));
t84 = t43 * t47;
t99 = -t46 * t84 + t77 * t50;
t51 = cos(qJ(2));
t48 = sin(qJ(1));
t68 = t48 * t77;
t87 = cos(qJ(1));
t24 = -t47 * t68 + t51 * t87;
t82 = t43 * t50;
t9 = -t24 * t46 + t48 * t82;
t42 = qJ(3) + pkin(10);
t39 = sin(t42);
t40 = cos(t42);
t98 = -m(4) * pkin(2) - t50 * mrSges(4,1) + t46 * mrSges(4,2) + t93 * t39 - t95 * t40 - mrSges(3,1);
t62 = t77 * t87;
t22 = t47 * t62 + t48 * t51;
t71 = t43 * t87;
t54 = t22 * t46 + t50 * t71;
t83 = t43 * t48;
t81 = t43 * t51;
t78 = t87 * pkin(1) + pkin(7) * t83;
t76 = t46 * t83;
t70 = -pkin(1) * t48 + pkin(7) * t71;
t4 = t22 * t40 - t39 * t71;
t3 = -t22 * t39 - t40 * t71;
t33 = t46 * t71;
t69 = -t22 * t50 + t33;
t23 = t47 * t87 + t51 * t68;
t38 = pkin(3) * t50 + pkin(2);
t44 = -qJ(4) - pkin(8);
t63 = pkin(3) * t76 - t23 * t44 + t24 * t38 + t78;
t21 = t47 * t48 - t51 * t62;
t16 = t39 * t77 + t40 * t84;
t10 = t24 * t50 + t76;
t8 = t24 * t40 + t39 * t83;
t7 = t24 * t39 - t40 * t83;
t2 = t23 * t45 + t49 * t8;
t1 = t23 * t49 - t45 * t8;
t5 = [(-t87 * mrSges(2,1) - m(3) * t78 - t24 * mrSges(3,1) - m(4) * (pkin(2) * t24 + t78) - t10 * mrSges(4,1) - t9 * mrSges(4,2) - m(5) * t63 - t8 * mrSges(5,1) - m(6) * (pkin(4) * t8 + t63) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t93 * t7 + (-mrSges(3,3) * t43 + mrSges(2,2)) * t48 + t59 * t23) * g(2) + (t48 * mrSges(2,1) + t87 * mrSges(2,2) - m(3) * t70 + t22 * mrSges(3,1) - mrSges(3,3) * t71 - m(4) * (-pkin(2) * t22 + t70) - t69 * mrSges(4,1) - t54 * mrSges(4,2) + t93 * t3 + t95 * t4 + t92 * t21 + t97 * (-pkin(3) * t33 - t21 * t44 + t22 * t38 - t70)) * g(1), (-t97 * (-t21 * t38 - t22 * t44) - t92 * t22 - t98 * t21) * g(2) + (-t97 * (-t23 * t38 - t24 * t44) - t92 * t24 - t98 * t23) * g(1) + (-t97 * t38 * t81 + (t98 * t51 + (t97 * t44 - t92) * t47) * t43) * g(3), (-t99 * mrSges(4,1) - (-t46 * t77 - t47 * t82) * mrSges(4,2) + t93 * t16 - t95 * (-t39 * t84 + t40 * t77)) * g(3) + (mrSges(4,1) * t54 - mrSges(4,2) * t69 - t95 * t3 + t4 * t93) * g(2) + (-t9 * mrSges(4,1) + t10 * mrSges(4,2) + t95 * t7 + t8 * t93) * g(1) + (-g(1) * t9 + g(2) * t54 - g(3) * t99) * t97 * pkin(3), t97 * (-g(1) * t23 - g(2) * t21 + g(3) * t81), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((t21 * t49 - t4 * t45) * mrSges(6,1) + (-t21 * t45 - t4 * t49) * mrSges(6,2)) - g(3) * ((-t16 * t45 - t49 * t81) * mrSges(6,1) + (-t16 * t49 + t45 * t81) * mrSges(6,2))];
taug = t5(:);

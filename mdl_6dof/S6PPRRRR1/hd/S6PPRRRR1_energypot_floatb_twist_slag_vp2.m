% Calculate potential energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:37
% EndTime: 2019-03-08 18:59:38
% DurationCPUTime: 0.93s
% Computational Cost: add. (547->123), mult. (1205->157), div. (0->0), fcn. (1486->16), ass. (0->59)
t52 = sin(pkin(7));
t56 = cos(pkin(7));
t57 = cos(pkin(6));
t53 = sin(pkin(6));
t54 = cos(pkin(13));
t84 = t53 * t54;
t31 = -t52 * t84 + t56 * t57;
t50 = sin(pkin(13));
t55 = cos(pkin(12));
t51 = sin(pkin(12));
t85 = t51 * t57;
t34 = -t50 * t55 - t54 * t85;
t82 = t53 * t56;
t24 = -t34 * t52 + t51 * t82;
t97 = -m(1) - m(2);
t96 = -m(6) - m(7);
t95 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t58 = sin(qJ(6));
t61 = cos(qJ(6));
t94 = -m(7) * pkin(5) - t61 * mrSges(7,1) + t58 * mrSges(7,2) - mrSges(6,1);
t93 = -m(5) * pkin(9) - t58 * mrSges(7,1) - t61 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t92 = cos(qJ(3));
t81 = t55 * t57;
t32 = -t50 * t51 + t54 * t81;
t23 = -t32 * t52 - t55 * t82;
t59 = sin(qJ(4));
t91 = t23 * t59;
t90 = t24 * t59;
t89 = t31 * t59;
t87 = t50 * t53;
t86 = t51 * t53;
t83 = t53 * t55;
t79 = qJ(2) * t53;
t76 = qJ(1) + r_base(3);
t75 = t52 * t92;
t74 = t56 * t92;
t73 = t55 * pkin(1) + t51 * t79 + r_base(1);
t72 = t57 * qJ(2) + t76;
t71 = t53 * t75;
t70 = t51 * pkin(1) - t55 * t79 + r_base(2);
t35 = -t50 * t85 + t54 * t55;
t69 = t35 * pkin(2) + t24 * pkin(8) + t73;
t68 = pkin(2) * t87 + t31 * pkin(8) + t72;
t33 = t50 * t81 + t51 * t54;
t65 = t33 * pkin(2) + pkin(8) * t23 + t70;
t63 = -pkin(10) - pkin(9);
t62 = cos(qJ(4));
t60 = sin(qJ(3));
t49 = qJ(4) + qJ(5);
t46 = cos(t49);
t45 = sin(t49);
t43 = pkin(4) * t62 + pkin(3);
t22 = t57 * t52 * t60 + (t54 * t56 * t60 + t50 * t92) * t53;
t21 = -t57 * t75 + t60 * t87 - t74 * t84;
t14 = t35 * t92 + (t34 * t56 + t52 * t86) * t60;
t13 = -t34 * t74 + t35 * t60 - t51 * t71;
t12 = t33 * t92 + (t32 * t56 - t52 * t83) * t60;
t11 = -t32 * t74 + t33 * t60 + t55 * t71;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t76 - mrSges(2,3) - m(3) * t72 - t57 * mrSges(3,3) - (t50 * mrSges(3,1) + t54 * mrSges(3,2)) * t53 - m(4) * t68 - t22 * mrSges(4,1) - t31 * mrSges(4,3) - m(5) * (pkin(3) * t22 + t68) - (t22 * t62 + t89) * mrSges(5,1) - (-t22 * t59 + t31 * t62) * mrSges(5,2) + t96 * (pkin(4) * t89 - t21 * t63 + t22 * t43 + t68) + t95 * (t22 * t45 - t31 * t46) + t94 * (t22 * t46 + t31 * t45) + t93 * t21) * g(3) + (-(t12 * t62 + t91) * mrSges(5,1) + mrSges(3,3) * t83 - m(3) * t70 - m(5) * (pkin(3) * t12 + t65) - m(4) * t65 - (-t12 * t59 + t23 * t62) * mrSges(5,2) - t55 * mrSges(2,2) - t51 * mrSges(2,1) - t32 * mrSges(3,2) - t33 * mrSges(3,1) - t23 * mrSges(4,3) - t12 * mrSges(4,1) - mrSges(1,2) + t97 * r_base(2) + t96 * (pkin(4) * t91 - t11 * t63 + t12 * t43 + t65) + t95 * (t12 * t45 - t23 * t46) + t94 * (t12 * t46 + t23 * t45) + t93 * t11) * g(2) + (-mrSges(3,3) * t86 - (t14 * t62 + t90) * mrSges(5,1) - m(3) * t73 - m(5) * (pkin(3) * t14 + t69) - m(4) * t69 - (-t14 * t59 + t24 * t62) * mrSges(5,2) - t55 * mrSges(2,1) + t51 * mrSges(2,2) - t34 * mrSges(3,2) - t35 * mrSges(3,1) - t24 * mrSges(4,3) - t14 * mrSges(4,1) - mrSges(1,1) + t97 * r_base(1) + t96 * (pkin(4) * t90 - t13 * t63 + t14 * t43 + t69) + t95 * (t14 * t45 - t24 * t46) + t94 * (t14 * t46 + t24 * t45) + t93 * t13) * g(1);
U  = t1;

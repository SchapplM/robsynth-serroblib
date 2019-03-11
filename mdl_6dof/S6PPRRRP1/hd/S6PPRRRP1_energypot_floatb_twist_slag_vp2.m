% Calculate potential energy for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:13
% EndTime: 2019-03-08 18:52:14
% DurationCPUTime: 0.85s
% Computational Cost: add. (572->123), mult. (1402->157), div. (0->0), fcn. (1753->14), ass. (0->61)
t46 = sin(pkin(7));
t50 = cos(pkin(7));
t51 = cos(pkin(6));
t47 = sin(pkin(6));
t48 = cos(pkin(12));
t84 = t47 * t48;
t67 = -t46 * t84 + t50 * t51;
t44 = sin(pkin(12));
t49 = cos(pkin(11));
t45 = sin(pkin(11));
t85 = t45 * t51;
t31 = -t44 * t49 - t48 * t85;
t82 = t47 * t50;
t68 = -t31 * t46 + t45 * t82;
t95 = -m(1) - m(2);
t94 = -mrSges(6,1) - mrSges(7,1);
t93 = -mrSges(6,2) - mrSges(7,2);
t53 = sin(qJ(5));
t92 = -m(7) * (pkin(5) * t53 + pkin(9)) + mrSges(4,2) - mrSges(5,3);
t91 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t90 = cos(qJ(3));
t89 = cos(qJ(4));
t87 = t44 * t47;
t86 = t45 * t47;
t83 = t47 * t49;
t81 = t49 * t51;
t79 = qJ(2) * t47;
t76 = qJ(1) + r_base(3);
t74 = t46 * t90;
t73 = t50 * t90;
t72 = t49 * pkin(1) + t45 * t79 + r_base(1);
t71 = t51 * qJ(2) + t76;
t70 = t47 * t74;
t29 = -t44 * t45 + t48 * t81;
t69 = -t29 * t46 - t49 * t82;
t66 = t45 * pkin(1) - t49 * t79 + r_base(2);
t32 = -t44 * t85 + t48 * t49;
t65 = t32 * pkin(2) + t68 * pkin(8) + t72;
t55 = sin(qJ(3));
t16 = t32 * t90 + (t31 * t50 + t46 * t86) * t55;
t64 = t16 * pkin(3) + t65;
t63 = pkin(2) * t87 + t67 * pkin(8) + t71;
t23 = t51 * t46 * t55 + (t48 * t50 * t55 + t44 * t90) * t47;
t62 = t23 * pkin(3) + t63;
t15 = -t31 * t73 + t32 * t55 - t45 * t70;
t61 = pkin(9) * t15 + t64;
t22 = -t51 * t74 + t55 * t87 - t73 * t84;
t60 = pkin(9) * t22 + t62;
t30 = t44 * t81 + t45 * t48;
t59 = t30 * pkin(2) + pkin(8) * t69 + t66;
t14 = t30 * t90 + (t29 * t50 - t46 * t83) * t55;
t58 = t14 * pkin(3) + t59;
t13 = -t29 * t73 + t30 * t55 + t49 * t70;
t57 = pkin(9) * t13 + t58;
t56 = cos(qJ(5));
t54 = sin(qJ(4));
t40 = pkin(5) * t56 + pkin(4);
t18 = t23 * t89 + t54 * t67;
t8 = t16 * t89 + t54 * t68;
t6 = t14 * t89 + t54 * t69;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t76 - mrSges(2,3) - m(3) * t71 - t51 * mrSges(3,3) - (t44 * mrSges(3,1) + t48 * mrSges(3,2)) * t47 - m(4) * t63 - t23 * mrSges(4,1) - t67 * mrSges(4,3) - m(5) * t60 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t60) - m(7) * (t18 * t40 + t62) + t93 * (-t18 * t53 + t22 * t56) + t92 * t22 + t94 * (t18 * t56 + t22 * t53) + t91 * (t23 * t54 - t67 * t89)) * g(3) + (mrSges(3,3) * t83 - m(7) * (t40 * t6 + t58) - m(3) * t66 - t69 * mrSges(4,3) - m(6) * (pkin(4) * t6 + t57) - m(5) * t57 - m(4) * t59 - t49 * mrSges(2,2) - t45 * mrSges(2,1) - t29 * mrSges(3,2) - t30 * mrSges(3,1) - t14 * mrSges(4,1) - t6 * mrSges(5,1) - mrSges(1,2) + t95 * r_base(2) + t94 * (t13 * t53 + t56 * t6) + t92 * t13 + t93 * (t13 * t56 - t6 * t53) + t91 * (t14 * t54 - t69 * t89)) * g(2) + (-m(7) * (t40 * t8 + t64) - m(3) * t72 - m(4) * t65 - t68 * mrSges(4,3) - m(6) * (pkin(4) * t8 + t61) - m(5) * t61 - t49 * mrSges(2,1) + t45 * mrSges(2,2) - t32 * mrSges(3,1) - t31 * mrSges(3,2) - t16 * mrSges(4,1) - t8 * mrSges(5,1) - mrSges(1,1) - mrSges(3,3) * t86 + t95 * r_base(1) + t94 * (t15 * t53 + t56 * t8) + t93 * (t15 * t56 - t8 * t53) + t92 * t15 + t91 * (t16 * t54 - t68 * t89)) * g(1);
U  = t1;

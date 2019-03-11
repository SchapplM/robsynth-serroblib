% Calculate potential energy for
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:38
% EndTime: 2019-03-10 02:32:39
% DurationCPUTime: 0.80s
% Computational Cost: add. (572->123), mult. (1402->153), div. (0->0), fcn. (1753->14), ass. (0->61)
t47 = cos(pkin(6));
t53 = sin(qJ(1));
t55 = cos(qJ(2));
t80 = t53 * t55;
t52 = sin(qJ(2));
t56 = cos(qJ(1));
t82 = t52 * t56;
t31 = -t47 * t80 - t82;
t44 = sin(pkin(7));
t46 = cos(pkin(7));
t45 = sin(pkin(6));
t86 = t45 * t53;
t68 = -t31 * t44 + t46 * t86;
t85 = t45 * t55;
t67 = -t44 * t85 + t46 * t47;
t95 = -m(1) - m(2);
t94 = -mrSges(6,1) - mrSges(7,1);
t93 = -mrSges(6,2) - mrSges(7,2);
t49 = sin(qJ(5));
t92 = -m(7) * (pkin(5) * t49 + pkin(11)) + mrSges(4,2) - mrSges(5,3);
t91 = -m(6) * pkin(12) + m(7) * (-qJ(6) - pkin(12)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t90 = cos(qJ(3));
t89 = cos(qJ(4));
t87 = t45 * t52;
t84 = t45 * t56;
t81 = t53 * t52;
t79 = t55 * t56;
t78 = pkin(8) + r_base(3);
t74 = t44 * t90;
t73 = t46 * t90;
t72 = pkin(9) * t47 + t78;
t71 = pkin(1) * t56 + pkin(9) * t86 + r_base(1);
t70 = t45 * t74;
t29 = t47 * t79 - t81;
t69 = -t29 * t44 - t46 * t84;
t66 = pkin(1) * t53 - pkin(9) * t84 + r_base(2);
t32 = -t47 * t81 + t79;
t65 = t32 * pkin(2) + pkin(10) * t68 + t71;
t64 = pkin(2) * t87 + pkin(10) * t67 + t72;
t51 = sin(qJ(3));
t18 = t32 * t90 + (t31 * t46 + t44 * t86) * t51;
t63 = pkin(3) * t18 + t65;
t23 = t47 * t44 * t51 + (t46 * t51 * t55 + t52 * t90) * t45;
t62 = pkin(3) * t23 + t64;
t17 = -t31 * t73 + t32 * t51 - t53 * t70;
t61 = pkin(11) * t17 + t63;
t22 = -t47 * t74 + t51 * t87 - t73 * t85;
t60 = pkin(11) * t22 + t62;
t30 = t47 * t82 + t80;
t59 = pkin(2) * t30 + pkin(10) * t69 + t66;
t16 = t30 * t90 + (t29 * t46 - t44 * t84) * t51;
t58 = pkin(3) * t16 + t59;
t15 = -t29 * t73 + t30 * t51 + t56 * t70;
t57 = t15 * pkin(11) + t58;
t54 = cos(qJ(5));
t50 = sin(qJ(4));
t40 = pkin(5) * t54 + pkin(4);
t14 = t23 * t89 + t50 * t67;
t10 = t18 * t89 + t50 * t68;
t8 = t16 * t89 + t50 * t69;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t78 - mrSges(2,3) - m(3) * t72 - t47 * mrSges(3,3) - (t52 * mrSges(3,1) + t55 * mrSges(3,2)) * t45 - m(4) * t64 - t23 * mrSges(4,1) - t67 * mrSges(4,3) - m(5) * t60 - t14 * mrSges(5,1) - m(6) * (pkin(4) * t14 + t60) - m(7) * (t14 * t40 + t62) + t94 * (t14 * t54 + t22 * t49) + t93 * (-t14 * t49 + t22 * t54) + t92 * t22 + t91 * (t23 * t50 - t67 * t89)) * g(3) + (mrSges(3,3) * t84 - m(7) * (t8 * t40 + t58) - t69 * mrSges(4,3) - m(3) * t66 - m(6) * (t8 * pkin(4) + t57) - m(5) * t57 - m(4) * t59 - t53 * mrSges(2,1) - t56 * mrSges(2,2) - t29 * mrSges(3,2) - t30 * mrSges(3,1) - t16 * mrSges(4,1) - t8 * mrSges(5,1) - mrSges(1,2) + t95 * r_base(2) + t94 * (t15 * t49 + t54 * t8) + t92 * t15 + t93 * (t15 * t54 - t49 * t8) + t91 * (t16 * t50 - t69 * t89)) * g(2) + (-m(7) * (t10 * t40 + t63) - m(3) * t71 - m(4) * t65 - t68 * mrSges(4,3) - m(6) * (pkin(4) * t10 + t61) - m(5) * t61 - mrSges(3,3) * t86 + t53 * mrSges(2,2) - t56 * mrSges(2,1) - t31 * mrSges(3,2) - t32 * mrSges(3,1) - t18 * mrSges(4,1) - t10 * mrSges(5,1) - mrSges(1,1) + t95 * r_base(1) + t94 * (t10 * t54 + t17 * t49) + t93 * (-t10 * t49 + t17 * t54) + t92 * t17 + t91 * (t18 * t50 - t68 * t89)) * g(1);
U  = t1;

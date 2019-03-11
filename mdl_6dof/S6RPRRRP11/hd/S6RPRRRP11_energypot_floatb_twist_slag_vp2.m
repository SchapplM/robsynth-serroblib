% Calculate potential energy for
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:21
% EndTime: 2019-03-09 06:34:22
% DurationCPUTime: 0.85s
% Computational Cost: add. (572->123), mult. (1402->155), div. (0->0), fcn. (1753->14), ass. (0->61)
t44 = sin(pkin(12));
t49 = cos(pkin(6));
t56 = cos(qJ(1));
t47 = cos(pkin(12));
t54 = sin(qJ(1));
t80 = t54 * t47;
t31 = -t44 * t56 - t49 * t80;
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t46 = sin(pkin(6));
t85 = t46 * t54;
t68 = -t31 * t45 + t48 * t85;
t86 = t46 * t47;
t67 = -t45 * t86 + t48 * t49;
t95 = -m(1) - m(2);
t94 = -mrSges(6,1) - mrSges(7,1);
t93 = -mrSges(6,2) - mrSges(7,2);
t51 = sin(qJ(5));
t92 = -m(7) * (pkin(5) * t51 + pkin(10)) + mrSges(4,2) - mrSges(5,3);
t91 = -m(6) * pkin(11) + m(7) * (-qJ(6) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t90 = cos(qJ(3));
t89 = cos(qJ(4));
t87 = t44 * t46;
t84 = t46 * t56;
t82 = t49 * t56;
t81 = t54 * t44;
t79 = qJ(2) * t46;
t78 = pkin(8) + r_base(3);
t74 = t45 * t90;
t73 = t48 * t90;
t72 = qJ(2) * t49 + t78;
t71 = pkin(1) * t56 + t54 * t79 + r_base(1);
t70 = t46 * t74;
t29 = t47 * t82 - t81;
t69 = -t29 * t45 - t48 * t84;
t66 = pkin(1) * t54 - t56 * t79 + r_base(2);
t32 = t47 * t56 - t49 * t81;
t65 = t32 * pkin(2) + pkin(9) * t68 + t71;
t64 = pkin(2) * t87 + pkin(9) * t67 + t72;
t53 = sin(qJ(3));
t18 = t32 * t90 + (t31 * t48 + t45 * t85) * t53;
t63 = pkin(3) * t18 + t65;
t23 = t49 * t45 * t53 + (t47 * t48 * t53 + t44 * t90) * t46;
t62 = t23 * pkin(3) + t64;
t17 = -t31 * t73 + t32 * t53 - t54 * t70;
t61 = pkin(10) * t17 + t63;
t22 = -t49 * t74 + t53 * t87 - t73 * t86;
t60 = pkin(10) * t22 + t62;
t30 = t44 * t82 + t80;
t59 = pkin(2) * t30 + pkin(9) * t69 + t66;
t16 = t30 * t90 + (t29 * t48 - t45 * t84) * t53;
t58 = pkin(3) * t16 + t59;
t15 = -t29 * t73 + t30 * t53 + t56 * t70;
t57 = t15 * pkin(10) + t58;
t55 = cos(qJ(5));
t52 = sin(qJ(4));
t40 = pkin(5) * t55 + pkin(4);
t14 = t23 * t89 + t52 * t67;
t10 = t18 * t89 + t52 * t68;
t8 = t16 * t89 + t52 * t69;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t78 - mrSges(2,3) - m(3) * t72 - t49 * mrSges(3,3) - (t44 * mrSges(3,1) + t47 * mrSges(3,2)) * t46 - m(4) * t64 - t23 * mrSges(4,1) - t67 * mrSges(4,3) - m(5) * t60 - t14 * mrSges(5,1) - m(6) * (pkin(4) * t14 + t60) - m(7) * (t14 * t40 + t62) + t94 * (t14 * t55 + t22 * t51) + t93 * (-t14 * t51 + t22 * t55) + t92 * t22 + t91 * (t23 * t52 - t67 * t89)) * g(3) + (mrSges(3,3) * t84 - m(7) * (t8 * t40 + t58) - m(3) * t66 - t69 * mrSges(4,3) - m(4) * t59 - m(6) * (t8 * pkin(4) + t57) - m(5) * t57 - t54 * mrSges(2,1) - t56 * mrSges(2,2) - t29 * mrSges(3,2) - t30 * mrSges(3,1) - t16 * mrSges(4,1) - t8 * mrSges(5,1) - mrSges(1,2) + t95 * r_base(2) + t94 * (t15 * t51 + t55 * t8) + t92 * t15 + t93 * (t15 * t55 - t51 * t8) + t91 * (t16 * t52 - t69 * t89)) * g(2) + (-m(7) * (t10 * t40 + t63) - m(3) * t71 - t68 * mrSges(4,3) - m(4) * t65 - m(6) * (pkin(4) * t10 + t61) - m(5) * t61 + t54 * mrSges(2,2) - t56 * mrSges(2,1) - t31 * mrSges(3,2) - t32 * mrSges(3,1) - t18 * mrSges(4,1) - t10 * mrSges(5,1) - mrSges(3,3) * t85 - mrSges(1,1) + t95 * r_base(1) + t94 * (t10 * t55 + t17 * t51) + t93 * (-t10 * t51 + t17 * t55) + t92 * t17 + t91 * (t18 * t52 - t68 * t89)) * g(1);
U  = t1;

% Calculate potential energy for
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:04
% EndTime: 2019-03-09 00:20:05
% DurationCPUTime: 0.84s
% Computational Cost: add. (572->123), mult. (1402->156), div. (0->0), fcn. (1753->14), ass. (0->60)
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t49 = cos(pkin(6));
t46 = sin(pkin(6));
t56 = cos(qJ(2));
t82 = t46 * t56;
t67 = -t45 * t82 + t49 * t48;
t44 = sin(pkin(12));
t47 = cos(pkin(12));
t54 = sin(qJ(2));
t79 = t49 * t56;
t31 = -t44 * t79 - t47 * t54;
t84 = t46 * t48;
t68 = -t31 * t45 + t44 * t84;
t94 = -m(1) - m(2);
t93 = -mrSges(6,1) - mrSges(7,1);
t92 = -mrSges(6,2) - mrSges(7,2);
t51 = sin(qJ(5));
t91 = -m(7) * (pkin(5) * t51 + pkin(10)) + mrSges(4,2) - mrSges(5,3);
t90 = -m(6) * pkin(11) + m(7) * (-qJ(6) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t89 = cos(qJ(3));
t88 = cos(qJ(4));
t86 = t44 * t46;
t85 = t46 * t47;
t83 = t46 * t54;
t80 = t49 * t54;
t76 = qJ(1) + r_base(3);
t74 = t45 * t89;
t73 = t48 * t89;
t72 = t47 * pkin(1) + pkin(8) * t86 + r_base(1);
t71 = t49 * pkin(8) + t76;
t70 = t46 * t74;
t29 = -t44 * t54 + t47 * t79;
t69 = -t29 * t45 - t47 * t84;
t66 = t44 * pkin(1) - pkin(8) * t85 + r_base(2);
t32 = -t44 * t80 + t47 * t56;
t65 = t32 * pkin(2) + t68 * pkin(9) + t72;
t53 = sin(qJ(3));
t16 = t32 * t89 + (t31 * t48 + t45 * t86) * t53;
t64 = t16 * pkin(3) + t65;
t63 = pkin(2) * t83 + t67 * pkin(9) + t71;
t23 = t49 * t45 * t53 + (t48 * t53 * t56 + t54 * t89) * t46;
t62 = t23 * pkin(3) + t63;
t15 = -t31 * t73 + t32 * t53 - t44 * t70;
t61 = pkin(10) * t15 + t64;
t22 = -t49 * t74 + t53 * t83 - t73 * t82;
t60 = t22 * pkin(10) + t62;
t30 = t44 * t56 + t47 * t80;
t59 = t30 * pkin(2) + pkin(9) * t69 + t66;
t14 = t30 * t89 + (t29 * t48 - t45 * t85) * t53;
t58 = t14 * pkin(3) + t59;
t13 = -t29 * t73 + t30 * t53 + t47 * t70;
t57 = pkin(10) * t13 + t58;
t55 = cos(qJ(5));
t52 = sin(qJ(4));
t40 = pkin(5) * t55 + pkin(4);
t18 = t23 * t88 + t52 * t67;
t10 = t16 * t88 + t52 * t68;
t8 = t14 * t88 + t52 * t69;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t76 - mrSges(2,3) - m(3) * t71 - t49 * mrSges(3,3) - (t54 * mrSges(3,1) + t56 * mrSges(3,2)) * t46 - m(4) * t63 - t23 * mrSges(4,1) - t67 * mrSges(4,3) - m(5) * t60 - t18 * mrSges(5,1) - m(6) * (t18 * pkin(4) + t60) - m(7) * (t18 * t40 + t62) + t93 * (t18 * t55 + t22 * t51) + t92 * (-t18 * t51 + t22 * t55) + t91 * t22 + t90 * (t23 * t52 - t67 * t88)) * g(3) + (mrSges(3,3) * t85 - m(7) * (t40 * t8 + t58) - t69 * mrSges(4,3) - m(3) * t66 - m(4) * t59 - m(6) * (pkin(4) * t8 + t57) - m(5) * t57 - t47 * mrSges(2,2) - t44 * mrSges(2,1) - t30 * mrSges(3,1) - t29 * mrSges(3,2) - t14 * mrSges(4,1) - t8 * mrSges(5,1) - mrSges(1,2) + t94 * r_base(2) + t93 * (t13 * t51 + t55 * t8) + t91 * t13 + t92 * (t13 * t55 - t51 * t8) + t90 * (t14 * t52 - t69 * t88)) * g(2) + (-m(7) * (t10 * t40 + t64) - m(3) * t72 - t68 * mrSges(4,3) - m(4) * t65 - m(6) * (pkin(4) * t10 + t61) - m(5) * t61 - t47 * mrSges(2,1) + t44 * mrSges(2,2) - t31 * mrSges(3,2) - t32 * mrSges(3,1) - t16 * mrSges(4,1) - t10 * mrSges(5,1) - mrSges(3,3) * t86 - mrSges(1,1) + t94 * r_base(1) + t93 * (t10 * t55 + t15 * t51) + t92 * (-t10 * t51 + t15 * t55) + t91 * t15 + t90 * (t16 * t52 - t68 * t88)) * g(1);
U  = t1;

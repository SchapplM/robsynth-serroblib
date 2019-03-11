% Calculate potential energy for
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:44
% EndTime: 2019-03-09 07:27:45
% DurationCPUTime: 0.88s
% Computational Cost: add. (547->123), mult. (1205->156), div. (0->0), fcn. (1486->16), ass. (0->58)
t50 = sin(pkin(13));
t53 = cos(pkin(13));
t62 = cos(qJ(1));
t55 = cos(pkin(6));
t59 = sin(qJ(1));
t81 = t55 * t59;
t34 = -t50 * t62 - t53 * t81;
t51 = sin(pkin(7));
t54 = cos(pkin(7));
t52 = sin(pkin(6));
t84 = t52 * t59;
t24 = -t34 * t51 + t54 * t84;
t85 = t52 * t53;
t31 = -t51 * t85 + t54 * t55;
t96 = -m(1) - m(2);
t95 = -m(6) - m(7);
t94 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t56 = sin(qJ(6));
t60 = cos(qJ(6));
t93 = -m(7) * pkin(5) - t60 * mrSges(7,1) + t56 * mrSges(7,2) - mrSges(6,1);
t92 = -m(5) * pkin(10) - t56 * mrSges(7,1) - t60 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t91 = cos(qJ(3));
t80 = t55 * t62;
t32 = -t50 * t59 + t53 * t80;
t83 = t52 * t62;
t23 = -t32 * t51 - t54 * t83;
t57 = sin(qJ(4));
t90 = t23 * t57;
t89 = t24 * t57;
t88 = t31 * t57;
t86 = t50 * t52;
t79 = qJ(2) * t52;
t78 = pkin(8) + r_base(3);
t75 = t51 * t91;
t74 = t54 * t91;
t73 = t55 * qJ(2) + t78;
t72 = t62 * pkin(1) + t59 * t79 + r_base(1);
t71 = t52 * t75;
t70 = t59 * pkin(1) - t62 * t79 + r_base(2);
t35 = -t50 * t81 + t53 * t62;
t69 = t35 * pkin(2) + t24 * pkin(9) + t72;
t68 = pkin(2) * t86 + t31 * pkin(9) + t73;
t33 = t50 * t80 + t53 * t59;
t65 = t33 * pkin(2) + pkin(9) * t23 + t70;
t63 = -pkin(11) - pkin(10);
t61 = cos(qJ(4));
t58 = sin(qJ(3));
t49 = qJ(4) + qJ(5);
t46 = cos(t49);
t45 = sin(t49);
t43 = pkin(4) * t61 + pkin(3);
t22 = t55 * t51 * t58 + (t53 * t54 * t58 + t50 * t91) * t52;
t21 = -t55 * t75 + t58 * t86 - t74 * t85;
t14 = t35 * t91 + (t34 * t54 + t51 * t84) * t58;
t13 = -t34 * t74 + t35 * t58 - t59 * t71;
t12 = t33 * t91 + (t32 * t54 - t51 * t83) * t58;
t11 = -t32 * t74 + t33 * t58 + t62 * t71;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t78 - mrSges(2,3) - m(3) * t73 - t55 * mrSges(3,3) - (t50 * mrSges(3,1) + t53 * mrSges(3,2)) * t52 - m(4) * t68 - t22 * mrSges(4,1) - t31 * mrSges(4,3) - m(5) * (pkin(3) * t22 + t68) - (t22 * t61 + t88) * mrSges(5,1) - (-t22 * t57 + t31 * t61) * mrSges(5,2) + t95 * (pkin(4) * t88 - t21 * t63 + t22 * t43 + t68) + t94 * (t22 * t45 - t31 * t46) + t93 * (t22 * t46 + t31 * t45) + t92 * t21) * g(3) + (-(t12 * t61 + t90) * mrSges(5,1) + mrSges(3,3) * t83 - m(3) * t70 - m(5) * (pkin(3) * t12 + t65) - m(4) * t65 - (-t12 * t57 + t23 * t61) * mrSges(5,2) - t62 * mrSges(2,2) - t59 * mrSges(2,1) - t32 * mrSges(3,2) - t33 * mrSges(3,1) - t23 * mrSges(4,3) - t12 * mrSges(4,1) - mrSges(1,2) + t96 * r_base(2) + t95 * (pkin(4) * t90 - t11 * t63 + t12 * t43 + t65) + t94 * (t12 * t45 - t23 * t46) + t93 * (t12 * t46 + t23 * t45) + t92 * t11) * g(2) + (-mrSges(3,3) * t84 - (t14 * t61 + t89) * mrSges(5,1) - m(3) * t72 - m(5) * (pkin(3) * t14 + t69) - m(4) * t69 - (-t14 * t57 + t24 * t61) * mrSges(5,2) - t62 * mrSges(2,1) + t59 * mrSges(2,2) - t35 * mrSges(3,1) - t34 * mrSges(3,2) - t24 * mrSges(4,3) - t14 * mrSges(4,1) - mrSges(1,1) + t96 * r_base(1) + t95 * (pkin(4) * t89 - t13 * t63 + t14 * t43 + t69) + t94 * (t14 * t45 - t24 * t46) + t93 * (t14 * t46 + t24 * t45) + t92 * t13) * g(1);
U  = t1;

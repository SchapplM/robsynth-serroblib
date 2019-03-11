% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:27
% EndTime: 2019-03-09 04:48:27
% DurationCPUTime: 0.15s
% Computational Cost: add. (148->68), mult. (203->79), div. (0->0), fcn. (203->8), ass. (0->45)
t79 = cos(qJ(4));
t66 = t79 * pkin(4) + pkin(3);
t75 = -qJ(5) - pkin(8);
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t106 = t66 * t77 + t75 * t80;
t105 = g(3) * pkin(6);
t104 = pkin(2) + pkin(6);
t81 = cos(qJ(1));
t103 = g(2) * t81;
t102 = g(3) * t80;
t74 = qJ(4) + pkin(9);
t67 = sin(t74);
t78 = sin(qJ(1));
t99 = t78 * t67;
t68 = cos(t74);
t98 = t78 * t68;
t76 = sin(qJ(4));
t97 = t78 * t76;
t96 = t78 * t79;
t95 = t81 * t67;
t94 = t81 * t68;
t93 = t81 * t76;
t92 = t81 * t79;
t70 = t78 * pkin(7);
t71 = t78 * pkin(1);
t91 = t70 + t71;
t90 = t81 * pkin(1) + t78 * qJ(2);
t89 = pkin(4) * t97 + t91;
t88 = t81 * pkin(7) + t90;
t87 = -t81 * qJ(2) + t71;
t86 = pkin(3) * t77 - pkin(8) * t80;
t59 = g(1) * t78 - t103;
t85 = t80 * t66 - t77 * t75 + t104;
t84 = pkin(4) * t93 + t106 * t78 + t88;
t53 = t77 * t99 - t94;
t55 = t77 * t95 + t98;
t83 = g(1) * t53 - g(2) * t55 + t67 * t102;
t82 = (-qJ(2) - t106) * t103;
t60 = g(1) * t81 + g(2) * t78;
t57 = -g(3) * t77 + t59 * t80;
t56 = -t77 * t94 + t99;
t54 = t77 * t98 + t95;
t52 = -g(1) * t54 - g(2) * t56 - t68 * t102;
t1 = [0, 0, 0, 0, 0, 0, -t60, t59, -g(3), -t105, 0, 0, 0, 0, 0, 0, -g(3), t60, -t59, -g(1) * t90 - g(2) * t87 - t105, 0, 0, 0, 0, 0, 0, -t59 * t77 - t102, -t57, -t60, -g(1) * t88 - g(2) * (t70 + t87) - g(3) * t104, 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t96 + t93) - g(2) * (-t77 * t92 + t97) - t79 * t102, -g(1) * (-t77 * t97 + t92) - g(2) * (t77 * t93 + t96) + t76 * t102, t57, -g(1) * (t86 * t78 + t88) - g(2) * t91 - g(3) * (t80 * pkin(3) + t77 * pkin(8) + t104) - (-qJ(2) - t86) * t103, 0, 0, 0, 0, 0, 0, t52, t83, t57, -g(1) * t84 - g(2) * t89 - g(3) * t85 - t82, 0, 0, 0, 0, 0, 0, t52, t57, -t83, -g(1) * (t54 * pkin(5) + t53 * qJ(6) + t84) - g(2) * (t56 * pkin(5) - t55 * qJ(6) + t89) - g(3) * ((pkin(5) * t68 + qJ(6) * t67) * t80 + t85) - t82;];
U_reg  = t1;

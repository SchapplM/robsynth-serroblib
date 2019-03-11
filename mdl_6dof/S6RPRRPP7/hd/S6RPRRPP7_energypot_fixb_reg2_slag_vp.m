% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:11
% EndTime: 2019-03-09 04:52:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (120->62), mult. (224->68), div. (0->0), fcn. (230->6), ass. (0->42)
t73 = sin(qJ(3));
t72 = sin(qJ(4));
t77 = cos(qJ(1));
t92 = t77 * t72;
t74 = sin(qJ(1));
t75 = cos(qJ(4));
t96 = t74 * t75;
t54 = t73 * t92 + t96;
t91 = t77 * t75;
t97 = t74 * t72;
t55 = -t73 * t91 + t97;
t103 = t55 * pkin(4) - t54 * qJ(5);
t102 = g(3) * pkin(6);
t101 = pkin(2) + pkin(6);
t100 = pkin(3) * t73;
t99 = g(2) * t77;
t76 = cos(qJ(3));
t98 = t72 * t76;
t95 = t74 * t76;
t94 = t75 * t76;
t93 = t76 * t77;
t90 = t77 * pkin(1) + t74 * qJ(2);
t88 = pkin(8) * t95;
t66 = t74 * pkin(7);
t67 = t74 * pkin(1);
t87 = pkin(8) * t93 + t66 + t67;
t86 = t77 * pkin(7) + t90;
t85 = t76 * pkin(3) + t73 * pkin(8) + t101;
t84 = -qJ(2) - t100;
t83 = -t77 * qJ(2) + t67;
t82 = t74 * t100 + t86;
t56 = g(1) * t74 - t99;
t81 = pkin(4) * t94 + qJ(5) * t98 + t85;
t52 = t73 * t97 - t91;
t53 = t73 * t96 + t92;
t80 = t53 * pkin(4) + t52 * qJ(5) + t82;
t79 = g(1) * t52 - g(2) * t54 + g(3) * t98;
t78 = t84 * t77 + t87;
t57 = g(1) * t77 + g(2) * t74;
t49 = g(1) * t95 - g(2) * t93 - g(3) * t73;
t48 = -g(1) * t53 - g(2) * t55 - g(3) * t94;
t1 = [0, 0, 0, 0, 0, 0, -t57, t56, -g(3), -t102, 0, 0, 0, 0, 0, 0, -g(3), t57, -t56, -g(1) * t90 - g(2) * t83 - t102, 0, 0, 0, 0, 0, 0, -g(3) * t76 - t56 * t73, -t49, -t57, -g(1) * t86 - g(2) * (t66 + t83) - g(3) * t101, 0, 0, 0, 0, 0, 0, t48, t79, t49, -g(1) * (t82 - t88) - g(2) * t78 - g(3) * t85, 0, 0, 0, 0, 0, 0, t48, t49, -t79, -g(1) * (t80 - t88) - g(2) * (t78 + t103) - g(3) * t81, 0, 0, 0, 0, 0, 0, t48, -t79, -t49, -g(1) * (t53 * pkin(5) + (-pkin(8) + qJ(6)) * t95 + t80) - g(2) * (t55 * pkin(5) + t103 + t87) - g(3) * (pkin(5) * t94 - t73 * qJ(6) + t81) - (-qJ(6) * t76 + t84) * t99;];
U_reg  = t1;

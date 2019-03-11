% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:26
% EndTime: 2019-03-09 03:46:26
% DurationCPUTime: 0.15s
% Computational Cost: add. (188->73), mult. (178->90), div. (0->0), fcn. (170->10), ass. (0->40)
t79 = qJ(1) + pkin(10);
t72 = sin(t79);
t86 = cos(qJ(3));
t106 = t72 * t86;
t83 = sin(qJ(3));
t98 = qJ(4) * t83;
t109 = pkin(3) * t106 + t72 * t98;
t81 = qJ(2) + pkin(6);
t108 = g(3) * t81;
t107 = g(3) * t86;
t73 = cos(t79);
t105 = t73 * t86;
t80 = qJ(5) + qJ(6);
t74 = sin(t80);
t104 = t74 * t83;
t75 = cos(t80);
t103 = t75 * t83;
t82 = sin(qJ(5));
t102 = t82 * t83;
t85 = cos(qJ(5));
t101 = t83 * t85;
t88 = -pkin(9) - pkin(8);
t100 = t86 * t88;
t84 = sin(qJ(1));
t99 = t84 * pkin(1) + t72 * pkin(2);
t97 = t72 * t102;
t96 = t83 * pkin(3) + t81;
t87 = cos(qJ(1));
t95 = t87 * pkin(1) + t73 * pkin(2) + t72 * pkin(7);
t94 = t99 + t109;
t93 = -t73 * pkin(7) + t99;
t92 = pkin(3) * t105 + t73 * t98 + t95;
t91 = g(1) * t73 + g(2) * t72;
t90 = -g(1) * t87 - g(2) * t84;
t89 = -t86 * qJ(4) + t96;
t71 = t85 * pkin(5) + pkin(4);
t62 = g(1) * t72 - g(2) * t73;
t61 = g(3) * t83 + t91 * t86;
t60 = t91 * t83 - t107;
t1 = [0, 0, 0, 0, 0, 0, t90, g(1) * t84 - g(2) * t87, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t91, t62, -g(3), t90 * pkin(1) - t108, 0, 0, 0, 0, 0, 0, -t61, t60, -t62, -g(1) * t95 - g(2) * t93 - t108, 0, 0, 0, 0, 0, 0, -t62, t61, -t60, -g(1) * t92 - g(2) * (t93 + t109) - g(3) * t89, 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t102 + t72 * t85) - g(2) * (-t73 * t85 + t97) + t82 * t107, -g(1) * (t73 * t101 - t72 * t82) - g(2) * (t72 * t101 + t73 * t82) + t85 * t107, -t61, -g(1) * (t72 * pkin(4) + pkin(8) * t105 + t92) - g(2) * (pkin(8) * t106 + (-pkin(4) - pkin(7)) * t73 + t94) - g(3) * (t83 * pkin(8) + t89) 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t104 + t72 * t75) - g(2) * (t72 * t104 - t73 * t75) + t74 * t107, -g(1) * (t73 * t103 - t72 * t74) - g(2) * (t72 * t103 + t73 * t74) + t75 * t107, -t61, -g(1) * (t72 * t71 + t92) - g(2) * (pkin(5) * t97 - t72 * t100 + t94) - g(3) * (-t83 * t88 + (-pkin(5) * t82 - qJ(4)) * t86 + t96) + (-g(1) * (pkin(5) * t102 - t100) - g(2) * (-pkin(7) - t71)) * t73;];
U_reg  = t1;

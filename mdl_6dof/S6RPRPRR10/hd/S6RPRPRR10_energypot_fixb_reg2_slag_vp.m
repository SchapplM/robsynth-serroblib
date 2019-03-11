% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:18
% EndTime: 2019-03-09 04:09:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->82), mult. (188->104), div. (0->0), fcn. (184->10), ass. (0->42)
t102 = g(3) * pkin(6);
t101 = pkin(2) + pkin(6);
t77 = sin(pkin(10));
t100 = pkin(4) * t77;
t83 = cos(qJ(1));
t99 = g(2) * t83;
t82 = cos(qJ(3));
t98 = g(3) * t82;
t78 = cos(pkin(10));
t65 = t78 * pkin(4) + pkin(3);
t80 = sin(qJ(3));
t81 = sin(qJ(1));
t97 = t80 * t81;
t96 = t80 * t83;
t76 = pkin(10) + qJ(5);
t68 = qJ(6) + t76;
t63 = sin(t68);
t95 = t81 * t63;
t64 = cos(t68);
t94 = t81 * t64;
t66 = sin(t76);
t93 = t81 * t66;
t67 = cos(t76);
t92 = t81 * t67;
t91 = t81 * t77;
t90 = t81 * t78;
t89 = t81 * t82;
t79 = -pkin(8) - qJ(4);
t71 = t81 * pkin(7);
t72 = t81 * pkin(1);
t88 = t71 + t72;
t87 = t83 * pkin(1) + t81 * qJ(2);
t86 = t83 * pkin(7) + t87;
t85 = -t83 * qJ(2) + t72;
t61 = g(1) * t81 - t99;
t84 = pkin(3) * t80 - qJ(4) * t82;
t75 = -pkin(9) + t79;
t62 = g(1) * t83 + g(2) * t81;
t60 = pkin(5) * t66 + t100;
t59 = pkin(5) * t67 + t65;
t58 = -g(3) * t80 + t61 * t82;
t1 = [0, 0, 0, 0, 0, 0, -t62, t61, -g(3), -t102, 0, 0, 0, 0, 0, 0, -g(3), t62, -t61, -g(1) * t87 - g(2) * t85 - t102, 0, 0, 0, 0, 0, 0, -t61 * t80 - t98, -t58, -t62, -g(1) * t86 - g(2) * (t71 + t85) - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t83 + t80 * t90) - g(2) * (-t78 * t96 + t91) - t78 * t98, -g(1) * (t78 * t83 - t80 * t91) - g(2) * (t77 * t96 + t90) + t77 * t98, t58, -g(1) * (t84 * t81 + t86) - g(2) * t88 - g(3) * (pkin(3) * t82 + qJ(4) * t80 + t101) - (-qJ(2) - t84) * t99, 0, 0, 0, 0, 0, 0, -g(1) * (t66 * t83 + t80 * t92) - g(2) * (-t67 * t96 + t93) - t67 * t98, -g(1) * (t67 * t83 - t80 * t93) - g(2) * (t66 * t96 + t92) + t66 * t98, t58, -g(1) * (t65 * t97 + t79 * t89 + t86) - g(2) * (pkin(4) * t91 + t88) - g(3) * (t65 * t82 - t79 * t80 + t101) + (-g(1) * t100 - g(2) * (-t65 * t80 - t79 * t82 - qJ(2))) * t83, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t83 + t80 * t94) - g(2) * (-t64 * t96 + t95) - t64 * t98, -g(1) * (t64 * t83 - t80 * t95) - g(2) * (t63 * t96 + t94) + t63 * t98, t58, -g(1) * (t59 * t97 + t75 * t89 + t86) - g(2) * (t81 * t60 + t88) - g(3) * (t59 * t82 - t75 * t80 + t101) + (-g(1) * t60 - g(2) * (-t59 * t80 - t75 * t82 - qJ(2))) * t83;];
U_reg  = t1;

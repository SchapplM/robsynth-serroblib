% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:14
% EndTime: 2019-03-09 11:41:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (201->63), mult. (174->78), div. (0->0), fcn. (166->10), ass. (0->37)
t90 = qJ(2) + pkin(10);
t85 = qJ(4) + t90;
t79 = sin(t85);
t80 = cos(t85);
t96 = cos(qJ(5));
t81 = t96 * pkin(5) + pkin(4);
t91 = -qJ(6) - pkin(9);
t113 = -t79 * t91 + t80 * t81;
t112 = g(3) * pkin(6);
t111 = g(3) * t79;
t94 = sin(qJ(2));
t110 = t94 * pkin(2) + pkin(6);
t97 = cos(qJ(2));
t82 = t97 * pkin(2) + pkin(1);
t93 = sin(qJ(5));
t95 = sin(qJ(1));
t107 = t95 * t93;
t106 = t95 * t96;
t98 = cos(qJ(1));
t105 = t98 * t93;
t104 = t98 * t96;
t92 = -pkin(7) - qJ(3);
t84 = cos(t90);
t74 = pkin(3) * t84 + t82;
t89 = -pkin(8) + t92;
t103 = t95 * t74 + t98 * t89;
t83 = sin(t90);
t102 = pkin(3) * t83 + t110;
t73 = t98 * t74;
t101 = -t95 * t89 + t73;
t100 = pkin(4) * t80 + pkin(9) * t79;
t99 = g(1) * t98 + g(2) * t95;
t75 = g(1) * t95 - g(2) * t98;
t71 = -g(3) * t80 + t99 * t79;
t70 = -g(1) * (t80 * t104 + t107) - g(2) * (t80 * t106 - t105) - t96 * t111;
t69 = -g(1) * (-t80 * t105 + t106) - g(2) * (-t80 * t107 - t104) + t93 * t111;
t1 = [0, 0, 0, 0, 0, 0, -t99, t75, -g(3), -t112, 0, 0, 0, 0, 0, 0, -g(3) * t94 - t99 * t97, -g(3) * t97 + t99 * t94, -t75, -g(1) * (t98 * pkin(1) + t95 * pkin(7)) - g(2) * (t95 * pkin(1) - t98 * pkin(7)) - t112, 0, 0, 0, 0, 0, 0, -g(3) * t83 - t99 * t84, -g(3) * t84 + t99 * t83, -t75, -g(1) * (t98 * t82 - t95 * t92) - g(2) * (t95 * t82 + t98 * t92) - g(3) * t110, 0, 0, 0, 0, 0, 0, -t99 * t80 - t111, t71, -t75, -g(1) * t101 - g(2) * t103 - g(3) * t102, 0, 0, 0, 0, 0, 0, t70, t69, -t71, -g(1) * (t100 * t98 + t101) - g(2) * (t100 * t95 + t103) - g(3) * (t79 * pkin(4) - t80 * pkin(9) + t102) 0, 0, 0, 0, 0, 0, t70, t69, -t71, -g(1) * (t113 * t98 + t73) - g(2) * (-pkin(5) * t105 + t103) - g(3) * (t79 * t81 + t80 * t91 + t102) + (-g(1) * (pkin(5) * t93 - t89) - g(2) * t113) * t95;];
U_reg  = t1;

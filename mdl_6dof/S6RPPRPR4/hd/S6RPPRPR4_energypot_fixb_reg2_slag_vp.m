% Calculate inertial parameters regressor of potential energy for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:07
% EndTime: 2019-03-09 01:47:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (152->57), mult. (217->72), div. (0->0), fcn. (247->10), ass. (0->35)
t101 = g(3) * pkin(6);
t76 = qJ(4) + pkin(10);
t69 = sin(t76);
t100 = g(3) * t69;
t78 = -qJ(3) + pkin(6);
t99 = g(3) * t78;
t98 = cos(qJ(1));
t97 = sin(qJ(1));
t70 = cos(t76);
t79 = sin(qJ(6));
t96 = t70 * t79;
t81 = cos(qJ(6));
t95 = t70 * t81;
t94 = t98 * pkin(1) + t97 * qJ(2);
t93 = cos(pkin(9));
t92 = sin(pkin(9));
t91 = t98 * pkin(2) + t94;
t80 = sin(qJ(4));
t90 = -t80 * pkin(4) + t78;
t61 = -t97 * t92 - t98 * t93;
t62 = t98 * t92 - t97 * t93;
t82 = cos(qJ(4));
t68 = t82 * pkin(4) + pkin(3);
t77 = -qJ(5) - pkin(7);
t89 = -t61 * t68 - t62 * t77 + t91;
t88 = -pkin(5) * t70 - pkin(8) * t69;
t87 = g(1) * t62 - g(2) * t61;
t86 = g(1) * t61 + g(2) * t62;
t85 = t97 * pkin(1) - t98 * qJ(2);
t84 = t97 * pkin(2) + t85;
t83 = t61 * t77 - t62 * t68 + t84;
t64 = -g(1) * t98 - g(2) * t97;
t63 = g(1) * t97 - g(2) * t98;
t55 = g(3) * t70 - t86 * t69;
t1 = [0, 0, 0, 0, 0, 0, t64, t63, -g(3), -t101, 0, 0, 0, 0, 0, 0, t64, -g(3), -t63, -g(1) * t94 - g(2) * t85 - t101, 0, 0, 0, 0, 0, 0, t86, t87, g(3), -g(1) * t91 - g(2) * t84 - t99, 0, 0, 0, 0, 0, 0, g(3) * t80 + t86 * t82, g(3) * t82 - t86 * t80, -t87, -g(1) * (-t61 * pkin(3) + t62 * pkin(7) + t91) - g(2) * (-t62 * pkin(3) - t61 * pkin(7) + t84) - t99, 0, 0, 0, 0, 0, 0, t86 * t70 + t100, t55, -t87, -g(1) * t89 - g(2) * t83 - g(3) * t90, 0, 0, 0, 0, 0, 0, -g(1) * (-t61 * t95 + t62 * t79) - g(2) * (-t61 * t79 - t62 * t95) + t81 * t100, -g(1) * (t61 * t96 + t62 * t81) - g(2) * (-t61 * t81 + t62 * t96) - t79 * t100, -t55, -g(1) * (t88 * t61 + t89) - g(2) * (t88 * t62 + t83) - g(3) * (-t69 * pkin(5) + t70 * pkin(8) + t90);];
U_reg  = t1;

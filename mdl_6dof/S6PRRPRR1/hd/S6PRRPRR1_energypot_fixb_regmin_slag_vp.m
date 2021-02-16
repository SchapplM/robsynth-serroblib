% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:21
% EndTime: 2021-01-16 03:28:22
% DurationCPUTime: 0.27s
% Computational Cost: add. (163->67), mult. (242->113), div. (0->0), fcn. (298->14), ass. (0->38)
t93 = sin(pkin(11));
t94 = sin(pkin(6));
t116 = t93 * t94;
t95 = cos(pkin(11));
t115 = t94 * t95;
t99 = sin(qJ(3));
t114 = t94 * t99;
t96 = cos(pkin(6));
t113 = t96 * t99;
t100 = sin(qJ(2));
t112 = t100 * t94;
t102 = cos(qJ(3));
t111 = t102 * t94;
t103 = cos(qJ(2));
t110 = t103 * t94;
t109 = t93 * t100;
t108 = t93 * t103;
t107 = t95 * t100;
t106 = t95 * t103;
t92 = qJ(3) + pkin(12);
t81 = -t96 * t106 + t109;
t83 = t96 * t108 + t107;
t104 = -g(1) * t83 - g(2) * t81 + g(3) * t110;
t101 = cos(qJ(6));
t98 = sin(qJ(6));
t97 = qJ(4) + pkin(8);
t91 = qJ(5) + t92;
t90 = cos(t92);
t89 = sin(t92);
t88 = t102 * pkin(3) + pkin(2);
t87 = cos(t91);
t86 = sin(t91);
t82 = t96 * t107 + t108;
t80 = t96 * t109 - t106;
t79 = t87 * t112 + t96 * t86;
t78 = t86 * t116 - t80 * t87;
t77 = -t86 * t115 + t82 * t87;
t1 = [-g(3) * qJ(1), 0, g(1) * t80 - g(2) * t82 - g(3) * t112, -t104, 0, 0, 0, 0, 0, -g(1) * (-t80 * t102 + t93 * t114) - g(2) * (t82 * t102 - t95 * t114) - g(3) * (t100 * t111 + t113), -g(1) * (t93 * t111 + t80 * t99) - g(2) * (-t95 * t111 - t82 * t99) - g(3) * (t96 * t102 - t99 * t112), -g(1) * (t89 * t116 - t80 * t90) - g(2) * (-t89 * t115 + t82 * t90) - g(3) * (t90 * t112 + t96 * t89), -g(1) * (t90 * t116 + t80 * t89) - g(2) * (-t90 * t115 - t82 * t89) - g(3) * (-t89 * t112 + t96 * t90), t104, -g(1) * (t95 * pkin(1) - t80 * t88 + t83 * t97) - g(2) * (t93 * pkin(1) + t81 * t97 + t82 * t88) - g(3) * (pkin(3) * t113 + t96 * pkin(7) + qJ(1)) + (-g(3) * (t100 * t88 - t103 * t97) + (-g(1) * t93 + g(2) * t95) * (pkin(3) * t99 + pkin(7))) * t94, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t77 - g(3) * t79, -g(1) * (t116 * t87 + t80 * t86) - g(2) * (-t115 * t87 - t82 * t86) - g(3) * (-t112 * t86 + t96 * t87), 0, 0, 0, 0, 0, -g(1) * (t78 * t101 + t83 * t98) - g(2) * (t77 * t101 + t81 * t98) - g(3) * (t79 * t101 - t110 * t98), -g(1) * (t83 * t101 - t78 * t98) - g(2) * (t81 * t101 - t77 * t98) - g(3) * (-t101 * t110 - t79 * t98);];
U_reg = t1;

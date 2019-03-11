% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:36:57
% EndTime: 2019-03-09 04:36:57
% DurationCPUTime: 0.14s
% Computational Cost: add. (209->61), mult. (222->71), div. (0->0), fcn. (228->8), ass. (0->38)
t82 = sin(qJ(3));
t84 = cos(qJ(4));
t100 = t82 * t84;
t81 = sin(qJ(4));
t102 = t81 * t82;
t107 = pkin(4) * t100 + qJ(5) * t102;
t85 = cos(qJ(3));
t106 = pkin(3) * t85;
t80 = qJ(2) + pkin(6);
t105 = g(3) * t80;
t79 = qJ(1) + pkin(9);
t73 = sin(t79);
t104 = t73 * t82;
t74 = cos(t79);
t103 = t74 * t82;
t101 = t81 * t85;
t99 = t84 * t85;
t98 = t82 * pkin(3) + t80;
t86 = cos(qJ(1));
t97 = t86 * pkin(1) + t74 * pkin(2) + t73 * pkin(7);
t83 = sin(qJ(1));
t96 = t83 * pkin(1) + t73 * pkin(2) - t74 * pkin(7);
t95 = pkin(8) * t103 + t74 * t106 + t97;
t94 = g(1) * t74 + g(2) * t73;
t93 = -g(1) * t86 - g(2) * t83;
t92 = -t85 * pkin(8) + t98;
t91 = pkin(8) * t104 + t73 * t106 + t96;
t57 = t73 * t101 + t74 * t84;
t59 = t74 * t101 - t73 * t84;
t90 = g(1) * t59 + g(2) * t57 + g(3) * t102;
t58 = t73 * t99 - t74 * t81;
t60 = t73 * t81 + t74 * t99;
t89 = g(1) * t60 + g(2) * t58 + g(3) * t100;
t88 = t60 * pkin(4) + t59 * qJ(5) + t95;
t87 = t58 * pkin(4) + t57 * qJ(5) + t91;
t61 = g(1) * t73 - g(2) * t74;
t54 = -g(3) * t85 + t94 * t82;
t1 = [0, 0, 0, 0, 0, 0, t93, g(1) * t83 - g(2) * t86, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t94, t61, -g(3), t93 * pkin(1) - t105, 0, 0, 0, 0, 0, 0, -g(3) * t82 - t94 * t85, t54, -t61, -g(1) * t97 - g(2) * t96 - t105, 0, 0, 0, 0, 0, 0, -t89, t90, -t54, -g(1) * t95 - g(2) * t91 - g(3) * t92, 0, 0, 0, 0, 0, 0, -t54, t89, -t90, -g(1) * t88 - g(2) * t87 - g(3) * (t92 + t107) 0, 0, 0, 0, 0, 0, -t54, -t90, -t89, -g(1) * (pkin(5) * t103 + t60 * qJ(6) + t88) - g(2) * (pkin(5) * t104 + t58 * qJ(6) + t87) - g(3) * (qJ(6) * t100 + (-pkin(5) - pkin(8)) * t85 + t98 + t107);];
U_reg  = t1;

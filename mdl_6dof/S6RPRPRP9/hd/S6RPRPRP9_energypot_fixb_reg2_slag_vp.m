% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:26
% EndTime: 2019-03-09 03:29:27
% DurationCPUTime: 0.15s
% Computational Cost: add. (148->68), mult. (203->79), div. (0->0), fcn. (203->8), ass. (0->45)
t78 = cos(pkin(9));
t68 = t78 * pkin(4) + pkin(3);
t79 = -pkin(8) - qJ(4);
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t108 = t68 * t80 + t79 * t82;
t107 = g(3) * pkin(6);
t106 = pkin(2) + pkin(6);
t83 = cos(qJ(1));
t105 = g(2) * t83;
t104 = g(3) * t82;
t76 = pkin(9) + qJ(5);
t69 = sin(t76);
t81 = sin(qJ(1));
t101 = t81 * t69;
t70 = cos(t76);
t100 = t81 * t70;
t77 = sin(pkin(9));
t99 = t81 * t77;
t98 = t81 * t78;
t97 = t83 * t69;
t96 = t83 * t70;
t95 = t83 * t77;
t94 = t83 * t78;
t72 = t81 * pkin(7);
t73 = t81 * pkin(1);
t93 = t72 + t73;
t92 = t83 * pkin(1) + t81 * qJ(2);
t91 = pkin(4) * t99 + t93;
t90 = t83 * pkin(7) + t92;
t89 = -t83 * qJ(2) + t73;
t62 = g(1) * t81 - t105;
t88 = pkin(3) * t80 - qJ(4) * t82;
t87 = t82 * t68 - t80 * t79 + t106;
t86 = pkin(4) * t95 + t108 * t81 + t90;
t55 = t80 * t101 - t96;
t57 = t80 * t97 + t100;
t85 = g(1) * t55 - g(2) * t57 + t69 * t104;
t84 = (-qJ(2) - t108) * t105;
t63 = g(1) * t83 + g(2) * t81;
t59 = -g(3) * t80 + t62 * t82;
t58 = -t80 * t96 + t101;
t56 = t80 * t100 + t97;
t54 = -g(1) * t56 - g(2) * t58 - t70 * t104;
t1 = [0, 0, 0, 0, 0, 0, -t63, t62, -g(3), -t107, 0, 0, 0, 0, 0, 0, -g(3), t63, -t62, -g(1) * t92 - g(2) * t89 - t107, 0, 0, 0, 0, 0, 0, -t62 * t80 - t104, -t59, -t63, -g(1) * t90 - g(2) * (t72 + t89) - g(3) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t98 + t95) - g(2) * (-t80 * t94 + t99) - t78 * t104, -g(1) * (-t80 * t99 + t94) - g(2) * (t80 * t95 + t98) + t77 * t104, t59, -g(1) * (t88 * t81 + t90) - g(2) * t93 - g(3) * (t82 * pkin(3) + t80 * qJ(4) + t106) - (-qJ(2) - t88) * t105, 0, 0, 0, 0, 0, 0, t54, t85, t59, -g(1) * t86 - g(2) * t91 - g(3) * t87 - t84, 0, 0, 0, 0, 0, 0, t54, t59, -t85, -g(1) * (t56 * pkin(5) + t55 * qJ(6) + t86) - g(2) * (t58 * pkin(5) - t57 * qJ(6) + t91) - g(3) * ((pkin(5) * t70 + qJ(6) * t69) * t82 + t87) - t84;];
U_reg  = t1;

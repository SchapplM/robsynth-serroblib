% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:31
% EndTime: 2019-03-09 06:32:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (148->68), mult. (203->79), div. (0->0), fcn. (203->8), ass. (0->45)
t79 = cos(qJ(4));
t67 = t79 * pkin(4) + pkin(3);
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t82 = -pkin(9) - pkin(8);
t107 = t67 * t77 + t80 * t82;
t106 = g(3) * pkin(6);
t105 = pkin(2) + pkin(6);
t81 = cos(qJ(1));
t104 = g(2) * t81;
t103 = g(3) * t80;
t75 = qJ(4) + qJ(5);
t68 = sin(t75);
t78 = sin(qJ(1));
t101 = t78 * t68;
t69 = cos(t75);
t100 = t78 * t69;
t76 = sin(qJ(4));
t99 = t78 * t76;
t98 = t78 * t79;
t96 = t81 * t68;
t95 = t81 * t69;
t94 = t81 * t76;
t93 = t81 * t79;
t71 = t78 * pkin(7);
t72 = t78 * pkin(1);
t92 = t71 + t72;
t91 = t81 * pkin(1) + t78 * qJ(2);
t90 = pkin(4) * t99 + t92;
t89 = t81 * pkin(7) + t91;
t88 = -t81 * qJ(2) + t72;
t87 = pkin(3) * t77 - pkin(8) * t80;
t60 = g(1) * t78 - t104;
t86 = t80 * t67 - t77 * t82 + t105;
t85 = pkin(4) * t94 + t107 * t78 + t89;
t54 = t77 * t101 - t95;
t56 = t77 * t96 + t100;
t84 = g(1) * t54 - g(2) * t56 + t68 * t103;
t83 = (-qJ(2) - t107) * t104;
t61 = g(1) * t81 + g(2) * t78;
t58 = -g(3) * t77 + t60 * t80;
t57 = -t77 * t95 + t101;
t55 = t77 * t100 + t96;
t53 = -g(1) * t55 - g(2) * t57 - t69 * t103;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t106, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t91 - g(2) * t88 - t106, 0, 0, 0, 0, 0, 0, -t60 * t77 - t103, -t58, -t61, -g(1) * t89 - g(2) * (t71 + t88) - g(3) * t105, 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t98 + t94) - g(2) * (-t77 * t93 + t99) - t79 * t103, -g(1) * (-t77 * t99 + t93) - g(2) * (t77 * t94 + t98) + t76 * t103, t58, -g(1) * (t87 * t78 + t89) - g(2) * t92 - g(3) * (t80 * pkin(3) + t77 * pkin(8) + t105) - (-qJ(2) - t87) * t104, 0, 0, 0, 0, 0, 0, t53, t84, t58, -g(1) * t85 - g(2) * t90 - g(3) * t86 - t83, 0, 0, 0, 0, 0, 0, t53, t58, -t84, -g(1) * (t55 * pkin(5) + t54 * qJ(6) + t85) - g(2) * (t57 * pkin(5) - t56 * qJ(6) + t90) - g(3) * ((pkin(5) * t69 + qJ(6) * t68) * t80 + t86) - t83;];
U_reg  = t1;

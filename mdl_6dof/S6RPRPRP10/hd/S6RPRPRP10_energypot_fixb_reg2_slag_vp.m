% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:31
% EndTime: 2019-03-09 03:32:31
% DurationCPUTime: 0.09s
% Computational Cost: add. (107->63), mult. (190->68), div. (0->0), fcn. (186->6), ass. (0->40)
t102 = g(3) * pkin(6);
t101 = pkin(2) + pkin(6);
t100 = -pkin(3) - pkin(8);
t82 = cos(qJ(1));
t99 = g(2) * t82;
t78 = sin(qJ(3));
t98 = g(3) * t78;
t79 = sin(qJ(1));
t97 = t78 * t79;
t81 = cos(qJ(3));
t96 = t79 * t81;
t77 = sin(qJ(5));
t95 = t82 * t77;
t80 = cos(qJ(5));
t94 = t82 * t80;
t93 = t82 * pkin(1) + t79 * qJ(2);
t92 = qJ(4) * t81;
t63 = t82 * t92;
t69 = t79 * pkin(7);
t71 = t79 * pkin(1);
t91 = t63 + t69 + t71;
t90 = t82 * pkin(7) + t93;
t89 = t81 * pkin(3) + t78 * qJ(4) + t101;
t88 = -t82 * qJ(2) + t71;
t87 = t69 + t88;
t60 = g(1) * t79 - t99;
t86 = g(3) * (t81 * pkin(8) + t89);
t85 = pkin(3) * t97 - t79 * t92 + t90;
t56 = t79 * t77 - t81 * t94;
t58 = t80 * t96 + t95;
t84 = -g(1) * t58 - g(2) * t56 + t80 * t98;
t83 = t82 * pkin(4) + pkin(8) * t97 + t85;
t70 = t79 * pkin(4);
t61 = g(1) * t82 + g(2) * t79;
t59 = -t77 * t96 + t94;
t57 = t79 * t80 + t81 * t95;
t55 = t60 * t81 - t98;
t54 = g(1) * t97 + g(3) * t81 - t78 * t99;
t53 = -g(1) * t59 - g(2) * t57 - t77 * t98;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t102, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t93 - g(2) * t88 - t102, 0, 0, 0, 0, 0, 0, -t54, -t55, -t61, -g(1) * t90 - g(2) * t87 - g(3) * t101, 0, 0, 0, 0, 0, 0, -t61, t54, t55, -g(1) * t85 - g(2) * ((-pkin(3) * t78 - qJ(2)) * t82 + t91) - g(3) * t89, 0, 0, 0, 0, 0, 0, t53, -t84, -t54, -g(1) * t83 - g(2) * (t70 + t91) - t86 - (t100 * t78 - qJ(2)) * t99, 0, 0, 0, 0, 0, 0, t53, -t54, t84, -g(1) * (t59 * pkin(5) + t58 * qJ(6) + t83) - g(2) * (t57 * pkin(5) + t56 * qJ(6) + t63 + t70 + t87) - t86 + (-g(3) * (pkin(5) * t77 - qJ(6) * t80) - t100 * t99) * t78;];
U_reg  = t1;

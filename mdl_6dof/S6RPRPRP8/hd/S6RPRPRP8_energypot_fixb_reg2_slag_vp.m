% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP8
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:03
% EndTime: 2019-03-09 03:26:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (148->61), mult. (181->69), div. (0->0), fcn. (177->8), ass. (0->41)
t103 = g(3) * pkin(6);
t102 = pkin(2) + pkin(6);
t78 = sin(qJ(3));
t101 = pkin(3) * t78;
t75 = qJ(3) + pkin(9);
t69 = sin(t75);
t100 = pkin(4) * t69;
t70 = cos(t75);
t99 = pkin(8) * t70;
t98 = g(3) * t70;
t77 = sin(qJ(5));
t79 = sin(qJ(1));
t97 = t79 * t77;
t80 = cos(qJ(5));
t96 = t79 * t80;
t82 = cos(qJ(1));
t95 = t82 * t77;
t94 = t82 * t80;
t93 = t82 * pkin(1) + t79 * qJ(2);
t81 = cos(qJ(3));
t92 = t81 * pkin(3) + t102;
t91 = t79 * t101 + t93;
t90 = -qJ(2) - t101;
t72 = t79 * pkin(1);
t76 = -qJ(4) - pkin(7);
t89 = -t79 * t76 + t72;
t88 = -t82 * qJ(2) + t72;
t87 = t70 * pkin(4) + t69 * pkin(8) + t92;
t86 = t82 * t99 + t89;
t61 = g(1) * t79 - g(2) * t82;
t85 = t91 + (t100 - t99) * t79;
t57 = t69 * t97 - t94;
t59 = t69 * t95 + t96;
t84 = g(1) * t57 - g(2) * t59 + t77 * t98;
t83 = (g(1) * t76 - g(2) * (t90 - t100)) * t82;
t62 = g(1) * t82 + g(2) * t79;
t60 = -t69 * t94 + t97;
t58 = t69 * t96 + t95;
t56 = -g(3) * t69 + t61 * t70;
t55 = -g(1) * t58 - g(2) * t60 - t80 * t98;
t1 = [0, 0, 0, 0, 0, 0, -t62, t61, -g(3), -t103, 0, 0, 0, 0, 0, 0, -g(3), t62, -t61, -g(1) * t93 - g(2) * t88 - t103, 0, 0, 0, 0, 0, 0, -g(3) * t81 - t61 * t78, g(3) * t78 - t61 * t81, -t62, -g(1) * (t82 * pkin(7) + t93) - g(2) * (t79 * pkin(7) + t88) - g(3) * t102, 0, 0, 0, 0, 0, 0, -t61 * t69 - t98, -t56, -t62, -g(1) * (-t82 * t76 + t91) - g(2) * (t90 * t82 + t89) - g(3) * t92, 0, 0, 0, 0, 0, 0, t55, t84, t56, -g(1) * t85 - g(2) * t86 - g(3) * t87 + t83, 0, 0, 0, 0, 0, 0, t55, t56, -t84, -g(1) * (t58 * pkin(5) + t57 * qJ(6) + t85) - g(2) * (t60 * pkin(5) - t59 * qJ(6) + t86) - g(3) * ((pkin(5) * t80 + qJ(6) * t77) * t70 + t87) + t83;];
U_reg  = t1;

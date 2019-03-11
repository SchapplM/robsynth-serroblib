% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:51
% EndTime: 2019-03-09 04:55:52
% DurationCPUTime: 0.15s
% Computational Cost: add. (120->62), mult. (224->67), div. (0->0), fcn. (230->6), ass. (0->41)
t71 = sin(qJ(3));
t70 = sin(qJ(4));
t75 = cos(qJ(1));
t91 = t75 * t70;
t72 = sin(qJ(1));
t73 = cos(qJ(4));
t94 = t72 * t73;
t53 = t71 * t91 + t94;
t90 = t75 * t73;
t95 = t72 * t70;
t54 = t71 * t90 - t95;
t101 = -t54 * pkin(4) - t53 * qJ(5);
t100 = g(3) * pkin(6);
t99 = pkin(2) + pkin(6);
t98 = pkin(3) * t71;
t97 = g(2) * t75;
t74 = cos(qJ(3));
t96 = t70 * t74;
t93 = t72 * t74;
t92 = t73 * t74;
t89 = t75 * pkin(1) + t72 * qJ(2);
t87 = pkin(8) * t93;
t65 = t72 * pkin(7);
t66 = t72 * pkin(1);
t86 = t75 * t74 * pkin(8) + t65 + t66;
t85 = t75 * pkin(7) + t89;
t84 = t74 * pkin(3) + t71 * pkin(8) + t99;
t83 = -qJ(2) - t98;
t82 = -t75 * qJ(2) + t66;
t81 = t72 * t98 + t85;
t55 = g(1) * t72 - t97;
t80 = pkin(4) * t92 + qJ(5) * t96 + t84;
t51 = t71 * t95 - t90;
t52 = t71 * t94 + t91;
t79 = t52 * pkin(4) + t51 * qJ(5) + t81;
t78 = g(1) * t51 - g(2) * t53 + g(3) * t96;
t77 = g(1) * t52 - g(2) * t54 + g(3) * t92;
t76 = t83 * t75 + t86;
t56 = g(1) * t75 + g(2) * t72;
t48 = -g(3) * t71 + t55 * t74;
t1 = [0, 0, 0, 0, 0, 0, -t56, t55, -g(3), -t100, 0, 0, 0, 0, 0, 0, -g(3), t56, -t55, -g(1) * t89 - g(2) * t82 - t100, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t55 * t71, -t48, -t56, -g(1) * t85 - g(2) * (t65 + t82) - g(3) * t99, 0, 0, 0, 0, 0, 0, -t77, t78, t48, -g(1) * (t81 - t87) - g(2) * t76 - g(3) * t84, 0, 0, 0, 0, 0, 0, t48, t77, -t78, -g(1) * (t79 - t87) - g(2) * (t76 + t101) - g(3) * t80, 0, 0, 0, 0, 0, 0, t48, -t78, -t77, -g(1) * (t52 * qJ(6) + (-pkin(5) - pkin(8)) * t93 + t79) - g(2) * (-t54 * qJ(6) + t101 + t86) - g(3) * (t71 * pkin(5) + qJ(6) * t92 + t80) - (pkin(5) * t74 + t83) * t97;];
U_reg  = t1;

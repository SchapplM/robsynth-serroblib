% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:46
% EndTime: 2019-03-09 02:26:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (160->66), mult. (261->90), div. (0->0), fcn. (305->10), ass. (0->40)
t101 = g(3) * pkin(6);
t72 = -qJ(3) + pkin(6);
t100 = g(3) * t72;
t74 = sin(qJ(4));
t99 = g(3) * t74;
t98 = cos(qJ(1));
t97 = sin(qJ(1));
t71 = qJ(5) + qJ(6);
t64 = sin(t71);
t76 = cos(qJ(4));
t96 = t64 * t76;
t65 = cos(t71);
t95 = t65 * t76;
t73 = sin(qJ(5));
t94 = t73 * t76;
t75 = cos(qJ(5));
t93 = t75 * t76;
t92 = t98 * pkin(1) + t97 * qJ(2);
t91 = cos(pkin(10));
t90 = sin(pkin(10));
t89 = t98 * pkin(2) + t92;
t88 = pkin(5) * t73 + pkin(7);
t56 = -t97 * t90 - t98 * t91;
t87 = -t56 * pkin(3) + t89;
t86 = -pkin(4) * t76 - pkin(8) * t74;
t57 = t98 * t90 - t97 * t91;
t85 = g(1) * t57 - g(2) * t56;
t84 = g(1) * t56 + g(2) * t57;
t83 = t97 * pkin(1) - t98 * qJ(2);
t63 = pkin(5) * t75 + pkin(4);
t77 = -pkin(9) - pkin(8);
t82 = t63 * t76 - t74 * t77;
t81 = t97 * pkin(2) + t83;
t80 = t57 * pkin(7) + t87;
t79 = -t57 * pkin(3) + t81;
t78 = -t56 * pkin(7) + t79;
t59 = -g(1) * t98 - g(2) * t97;
t58 = g(1) * t97 - g(2) * t98;
t53 = g(3) * t76 - t84 * t74;
t1 = [0, 0, 0, 0, 0, 0, t59, t58, -g(3), -t101, 0, 0, 0, 0, 0, 0, t59, -g(3), -t58, -g(1) * t92 - g(2) * t83 - t101, 0, 0, 0, 0, 0, 0, t84, t85, g(3), -g(1) * t89 - g(2) * t81 - t100, 0, 0, 0, 0, 0, 0, t84 * t76 + t99, t53, -t85, -g(1) * t80 - g(2) * t78 - t100, 0, 0, 0, 0, 0, 0, -g(1) * (-t56 * t93 + t57 * t73) - g(2) * (-t56 * t73 - t57 * t93) + t75 * t99, -g(1) * (t56 * t94 + t57 * t75) - g(2) * (-t56 * t75 + t57 * t94) - t73 * t99, -t53, -g(1) * (t86 * t56 + t80) - g(2) * (t86 * t57 + t78) - g(3) * (-pkin(4) * t74 + pkin(8) * t76 + t72) 0, 0, 0, 0, 0, 0, -g(1) * (-t56 * t95 + t57 * t64) - g(2) * (-t56 * t64 - t57 * t95) + t65 * t99, -g(1) * (t56 * t96 + t57 * t65) - g(2) * (-t56 * t65 + t57 * t96) - t64 * t99, -t53, -g(1) * t87 - g(2) * t79 - g(3) * (-t74 * t63 - t76 * t77 + t72) + (-g(1) * t88 + g(2) * t82) * t57 + (g(1) * t82 + g(2) * t88) * t56;];
U_reg  = t1;

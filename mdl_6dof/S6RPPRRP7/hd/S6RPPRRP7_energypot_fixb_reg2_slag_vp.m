% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:53
% EndTime: 2019-03-09 02:13:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (140->63), mult. (168->73), div. (0->0), fcn. (160->8), ass. (0->36)
t66 = pkin(9) + qJ(4);
t60 = sin(t66);
t61 = cos(t66);
t93 = pkin(4) * t60 - pkin(8) * t61;
t92 = g(3) * pkin(6);
t91 = pkin(2) + pkin(6);
t67 = sin(pkin(9));
t90 = pkin(3) * t67;
t87 = g(3) * t61;
t71 = sin(qJ(5));
t72 = sin(qJ(1));
t86 = t72 * t71;
t73 = cos(qJ(5));
t85 = t72 * t73;
t74 = cos(qJ(1));
t84 = t74 * t71;
t83 = t74 * t73;
t82 = t74 * pkin(1) + t72 * qJ(2);
t68 = cos(pkin(9));
t81 = t68 * pkin(3) + t91;
t80 = t72 * t90 + t82;
t70 = -pkin(7) - qJ(3);
t79 = pkin(5) * t71 - t70;
t78 = -qJ(2) - t90;
t64 = t72 * pkin(1);
t77 = -t72 * t70 + t64;
t76 = -t74 * qJ(2) + t64;
t56 = g(1) * t72 - g(2) * t74;
t59 = t73 * pkin(5) + pkin(4);
t69 = -qJ(6) - pkin(8);
t75 = t59 * t60 + t61 * t69;
t57 = g(1) * t74 + g(2) * t72;
t55 = -g(3) * t60 + t56 * t61;
t54 = -g(1) * (t60 * t85 + t84) - g(2) * (-t60 * t83 + t86) - t73 * t87;
t53 = -g(1) * (-t60 * t86 + t83) - g(2) * (t60 * t84 + t85) + t71 * t87;
t1 = [0, 0, 0, 0, 0, 0, -t57, t56, -g(3), -t92, 0, 0, 0, 0, 0, 0, -g(3), t57, -t56, -g(1) * t82 - g(2) * t76 - t92, 0, 0, 0, 0, 0, 0, -g(3) * t68 - t56 * t67, g(3) * t67 - t56 * t68, -t57, -g(1) * (t74 * qJ(3) + t82) - g(2) * (t72 * qJ(3) + t76) - g(3) * t91, 0, 0, 0, 0, 0, 0, -t56 * t60 - t87, -t55, -t57, -g(1) * (-t74 * t70 + t80) - g(2) * (t78 * t74 + t77) - g(3) * t81, 0, 0, 0, 0, 0, 0, t54, t53, t55, -g(1) * (t93 * t72 + t80) - g(2) * t77 - g(3) * (t61 * pkin(4) + t60 * pkin(8) + t81) + (g(1) * t70 - g(2) * (t78 - t93)) * t74, 0, 0, 0, 0, 0, 0, t54, t53, t55, -g(1) * t80 - g(2) * t64 - g(3) * (t61 * t59 - t60 * t69 + t81) + (-g(1) * t75 - g(2) * t79) * t72 + (-g(1) * t79 - g(2) * (-t75 + t78)) * t74;];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:57
% EndTime: 2019-03-09 01:44:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (167->58), mult. (129->66), div. (0->0), fcn. (117->10), ass. (0->35)
t66 = qJ(4) + pkin(10);
t59 = sin(t66);
t61 = cos(t66);
t94 = pkin(5) * t59 - pkin(8) * t61;
t71 = sin(qJ(4));
t93 = pkin(4) * t71;
t90 = g(3) * t61;
t69 = qJ(2) + pkin(6);
t89 = g(3) * t69;
t67 = qJ(1) + pkin(9);
t60 = sin(t67);
t70 = sin(qJ(6));
t88 = t60 * t70;
t73 = cos(qJ(6));
t87 = t60 * t73;
t62 = cos(t67);
t86 = t62 * t70;
t85 = t62 * t73;
t72 = sin(qJ(1));
t84 = t72 * pkin(1) + t60 * pkin(2);
t83 = pkin(3) + t69;
t75 = cos(qJ(1));
t82 = t75 * pkin(1) + t62 * pkin(2) + t60 * qJ(3);
t81 = -qJ(3) - t93;
t74 = cos(qJ(4));
t80 = t74 * pkin(4) + t83;
t79 = t60 * t93 + t82;
t68 = -qJ(5) - pkin(7);
t78 = -t60 * t68 + t84;
t77 = -t62 * qJ(3) + t84;
t53 = g(1) * t60 - g(2) * t62;
t76 = -g(1) * t75 - g(2) * t72;
t54 = g(1) * t62 + g(2) * t60;
t52 = -g(3) * t59 + t53 * t61;
t1 = [0, 0, 0, 0, 0, 0, t76, g(1) * t72 - g(2) * t75, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t54, t53, -g(3), t76 * pkin(1) - t89, 0, 0, 0, 0, 0, 0, -g(3), t54, -t53, -g(1) * t82 - g(2) * t77 - t89, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t53 * t71, g(3) * t71 - t53 * t74, -t54, -g(1) * (t62 * pkin(7) + t82) - g(2) * (t60 * pkin(7) + t77) - g(3) * t83, 0, 0, 0, 0, 0, 0, -t53 * t59 - t90, -t52, -t54, -g(1) * (-t62 * t68 + t79) - g(2) * (t81 * t62 + t78) - g(3) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (t59 * t87 + t86) - g(2) * (-t59 * t85 + t88) - t73 * t90, -g(1) * (-t59 * t88 + t85) - g(2) * (t59 * t86 + t87) + t70 * t90, t52, -g(1) * (t94 * t60 + t79) - g(2) * t78 - g(3) * (t61 * pkin(5) + t59 * pkin(8) + t80) + (g(1) * t68 - g(2) * (t81 - t94)) * t62;];
U_reg  = t1;

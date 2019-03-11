% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:28
% EndTime: 2019-03-09 02:31:28
% DurationCPUTime: 0.11s
% Computational Cost: add. (106->66), mult. (151->76), div. (0->0), fcn. (143->8), ass. (0->40)
t90 = g(3) * pkin(6);
t89 = pkin(2) + pkin(6);
t62 = sin(qJ(5));
t88 = pkin(5) * t62;
t66 = cos(qJ(4));
t87 = g(3) * t66;
t61 = qJ(5) + qJ(6);
t52 = sin(t61);
t64 = sin(qJ(1));
t86 = t64 * t52;
t53 = cos(t61);
t85 = t64 * t53;
t84 = t64 * t62;
t65 = cos(qJ(5));
t83 = t64 * t65;
t67 = cos(qJ(1));
t82 = t67 * t52;
t81 = t67 * t53;
t80 = t67 * t62;
t79 = t67 * t65;
t78 = t67 * pkin(1) + t64 * qJ(2);
t77 = pkin(3) + t89;
t76 = t67 * qJ(3) + t78;
t57 = t64 * pkin(1);
t75 = -t67 * qJ(2) + t57;
t74 = g(1) * t76;
t54 = t64 * qJ(3);
t73 = t54 + t75;
t63 = sin(qJ(4));
t72 = pkin(4) * t63 - pkin(8) * t66;
t50 = g(1) * t67 + g(2) * t64;
t51 = t65 * pkin(5) + pkin(4);
t68 = -pkin(9) - pkin(8);
t71 = t51 * t63 + t66 * t68;
t70 = -t64 * pkin(7) + t76;
t58 = t67 * pkin(7);
t69 = t58 + t73;
t49 = g(1) * t64 - g(2) * t67;
t48 = -g(3) * t63 + t50 * t66;
t1 = [0, 0, 0, 0, 0, 0, -t50, t49, -g(3), -t90, 0, 0, 0, 0, 0, 0, -g(3), t50, -t49, -g(1) * t78 - g(2) * t75 - t90, 0, 0, 0, 0, 0, 0, -g(3), -t49, -t50, -g(2) * t73 - g(3) * t89 - t74, 0, 0, 0, 0, 0, 0, -t50 * t63 - t87, -t48, t49, -g(1) * t70 - g(2) * t69 - g(3) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t79 - t84) - g(2) * (t63 * t83 + t80) - t65 * t87, -g(1) * (-t63 * t80 - t83) - g(2) * (-t63 * t84 + t79) + t62 * t87, t48, -g(1) * (t72 * t67 + t70) - g(2) * (t72 * t64 + t69) - g(3) * (t66 * pkin(4) + t63 * pkin(8) + t77) 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t81 - t86) - g(2) * (t63 * t85 + t82) - t53 * t87, -g(1) * (-t63 * t82 - t85) - g(2) * (-t63 * t86 + t81) + t52 * t87, t48, -t74 - g(2) * (t54 + t57 + t58) - g(3) * (t66 * t51 - t63 * t68 + t77) + (-g(1) * t71 - g(2) * (-qJ(2) + t88)) * t67 + (-g(1) * (-pkin(7) - t88) - g(2) * t71) * t64;];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:32
% EndTime: 2019-03-09 01:51:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (88->59), mult. (138->62), div. (0->0), fcn. (126->6), ass. (0->33)
t84 = g(3) * pkin(6);
t83 = pkin(2) + pkin(6);
t63 = sin(qJ(4));
t82 = pkin(4) * t63;
t81 = g(3) * t63;
t64 = sin(qJ(1));
t66 = cos(qJ(4));
t80 = t64 * t66;
t62 = sin(qJ(6));
t67 = cos(qJ(1));
t79 = t67 * t62;
t65 = cos(qJ(6));
t78 = t67 * t65;
t77 = t67 * pkin(1) + t64 * qJ(2);
t76 = qJ(5) * t66;
t75 = pkin(3) + t83;
t74 = t67 * qJ(3) + t77;
t57 = t64 * pkin(1);
t73 = -t67 * qJ(2) + t57;
t72 = t66 * pkin(4) + t63 * qJ(5) + t75;
t54 = t64 * qJ(3);
t71 = t54 + t73;
t50 = g(1) * t67 + g(2) * t64;
t70 = pkin(8) * t63 - t76;
t69 = -t64 * pkin(7) + t74;
t59 = t67 * pkin(7);
t68 = t59 + t71;
t52 = t67 * t82;
t51 = t64 * t82;
t49 = g(1) * t64 - g(2) * t67;
t48 = t50 * t66 - t81;
t47 = g(3) * t66 + t50 * t63;
t1 = [0, 0, 0, 0, 0, 0, -t50, t49, -g(3), -t84, 0, 0, 0, 0, 0, 0, -g(3), t50, -t49, -g(1) * t77 - g(2) * t73 - t84, 0, 0, 0, 0, 0, 0, -g(3), -t49, -t50, -g(1) * t74 - g(2) * t71 - g(3) * t83, 0, 0, 0, 0, 0, 0, -t47, -t48, t49, -g(1) * t69 - g(2) * t68 - g(3) * t75, 0, 0, 0, 0, 0, 0, t49, t47, t48, -g(1) * (-t67 * t76 + t52 + t69) - g(2) * (-t64 * t76 + t51 + t68) - g(3) * t72, 0, 0, 0, 0, 0, 0, -g(1) * (-t64 * t65 - t66 * t79) - g(2) * (-t62 * t80 + t78) - t62 * t81, -g(1) * (t64 * t62 - t66 * t78) - g(2) * (-t65 * t80 - t79) - t65 * t81, -t47, -g(1) * (t52 + t74) - g(2) * (t51 + t54 + t57 + t59) - g(3) * (t66 * pkin(8) + t72) + (-g(1) * t70 - g(2) * (pkin(5) - qJ(2))) * t67 + (-g(1) * (-pkin(5) - pkin(7)) - g(2) * t70) * t64;];
U_reg  = t1;

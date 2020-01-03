% Calculate inertial parameters regressor of potential energy for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:10
% EndTime: 2019-12-31 18:22:10
% DurationCPUTime: 0.11s
% Computational Cost: add. (138->59), mult. (130->79), div. (0->0), fcn. (125->10), ass. (0->32)
t65 = qJ(2) + pkin(5);
t84 = g(3) * t65;
t67 = sin(qJ(3));
t83 = g(3) * t67;
t62 = qJ(1) + pkin(8);
t55 = sin(t62);
t63 = sin(pkin(9));
t82 = t55 * t63;
t69 = cos(qJ(3));
t81 = t55 * t69;
t57 = cos(t62);
t80 = t57 * t69;
t79 = t63 * t69;
t64 = cos(pkin(9));
t78 = t64 * t69;
t66 = -pkin(7) - qJ(4);
t77 = t66 * t67;
t68 = sin(qJ(1));
t76 = t68 * pkin(1) + t55 * pkin(2);
t70 = cos(qJ(1));
t75 = t70 * pkin(1) + t57 * pkin(2) + t55 * pkin(6);
t74 = -t57 * pkin(6) + t76;
t73 = g(1) * t57 + g(2) * t55;
t72 = -g(1) * t70 - g(2) * t68;
t71 = pkin(3) * t69 + qJ(4) * t67;
t61 = pkin(9) + qJ(5);
t56 = cos(t61);
t54 = sin(t61);
t53 = t64 * pkin(4) + pkin(3);
t49 = g(1) * t55 - g(2) * t57;
t48 = -g(3) * t69 + t73 * t67;
t1 = [0, 0, 0, 0, 0, 0, t72, g(1) * t68 - g(2) * t70, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t73, t49, -g(3), t72 * pkin(1) - t84, 0, 0, 0, 0, 0, 0, -t73 * t69 - t83, t48, -t49, -g(1) * t75 - g(2) * t74 - t84, 0, 0, 0, 0, 0, 0, -g(1) * (t57 * t78 + t82) - g(2) * (t55 * t78 - t57 * t63) - t64 * t83, -g(1) * (t55 * t64 - t57 * t79) - g(2) * (-t55 * t79 - t57 * t64) + t63 * t83, -t48, -g(1) * (t71 * t57 + t75) - g(2) * (t71 * t55 + t74) - g(3) * (t67 * pkin(3) - t69 * qJ(4) + t65), 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t54 + t56 * t80) - g(2) * (-t57 * t54 + t56 * t81) - t56 * t83, -g(1) * (-t54 * t80 + t55 * t56) - g(2) * (-t54 * t81 - t57 * t56) + t54 * t83, -t48, -g(1) * (pkin(4) * t82 + t75) - g(2) * (t53 * t81 - t55 * t77 + t76) - g(3) * (t67 * t53 + t69 * t66 + t65) + (-g(1) * (t53 * t69 - t77) - g(2) * (-pkin(4) * t63 - pkin(6))) * t57;];
U_reg = t1;

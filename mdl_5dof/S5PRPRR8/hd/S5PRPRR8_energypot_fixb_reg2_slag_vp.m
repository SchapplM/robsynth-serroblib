% Calculate inertial parameters regressor of potential energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:25
% EndTime: 2019-12-05 16:04:25
% DurationCPUTime: 0.17s
% Computational Cost: add. (155->68), mult. (347->100), div. (0->0), fcn. (405->10), ass. (0->44)
t68 = sin(pkin(9));
t69 = sin(pkin(5));
t98 = t68 * t69;
t70 = cos(pkin(9));
t97 = t69 * t70;
t73 = sin(qJ(4));
t96 = t69 * t73;
t74 = sin(qJ(2));
t95 = t69 * t74;
t76 = cos(qJ(4));
t94 = t69 * t76;
t77 = cos(qJ(2));
t93 = t69 * t77;
t71 = cos(pkin(5));
t92 = t71 * t74;
t91 = t71 * t77;
t90 = t70 * pkin(1) + pkin(6) * t98;
t89 = t71 * pkin(6) + qJ(1);
t88 = pkin(6) * t97;
t52 = t68 * t74 - t70 * t91;
t53 = t68 * t77 + t70 * t92;
t64 = t68 * pkin(1);
t87 = t53 * pkin(2) + t52 * qJ(3) + t64;
t86 = g(1) * t68 - g(2) * t70;
t54 = t68 * t91 + t70 * t74;
t55 = -t68 * t92 + t70 * t77;
t85 = t55 * pkin(2) + t54 * qJ(3) + t90;
t84 = pkin(2) * t95 - qJ(3) * t93 + t89;
t44 = -t54 * t76 + t68 * t96;
t46 = t52 * t76 + t70 * t96;
t56 = t71 * t73 + t76 * t93;
t83 = g(1) * t44 - g(2) * t46 + g(3) * t56;
t82 = -g(1) * t54 - g(2) * t52 + g(3) * t93;
t81 = g(1) * t55 + g(2) * t53 + g(3) * t95;
t80 = t71 * pkin(3) + pkin(7) * t95 + t84;
t79 = pkin(3) * t98 + t55 * pkin(7) + t85;
t78 = t53 * pkin(7) + (-pkin(3) - pkin(6)) * t97 + t87;
t75 = cos(qJ(5));
t72 = sin(qJ(5));
t57 = t71 * t76 - t73 * t93;
t48 = -g(3) * t71 - t86 * t69;
t47 = t52 * t73 - t70 * t94;
t45 = t54 * t73 + t68 * t94;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68, t86, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t81, -t82, t48, -g(1) * t90 - g(2) * (t64 - t88) - g(3) * t89, 0, 0, 0, 0, 0, 0, t48, t81, t82, -g(1) * t85 - g(2) * (t87 - t88) - g(3) * t84, 0, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t47 - g(3) * t57, t83, -t81, -g(1) * t79 - g(2) * t78 - g(3) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (t45 * t75 + t55 * t72) - g(2) * (t47 * t75 + t53 * t72) - g(3) * (t57 * t75 + t72 * t95), -g(1) * (-t45 * t72 + t55 * t75) - g(2) * (-t47 * t72 + t53 * t75) - g(3) * (-t57 * t72 + t75 * t95), -t83, -g(1) * (t45 * pkin(4) + t44 * pkin(8) + t79) - g(2) * (t47 * pkin(4) - t46 * pkin(8) + t78) - g(3) * (t57 * pkin(4) + t56 * pkin(8) + t80);];
U_reg = t1;

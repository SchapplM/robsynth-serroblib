% Calculate inertial parameters regressor of potential energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:34
% EndTime: 2019-12-05 15:51:34
% DurationCPUTime: 0.26s
% Computational Cost: add. (233->75), mult. (538->124), div. (0->0), fcn. (668->12), ass. (0->48)
t76 = sin(pkin(10));
t79 = cos(pkin(10));
t84 = sin(qJ(2));
t87 = cos(qJ(2));
t66 = -t84 * t76 + t87 * t79;
t78 = sin(pkin(5));
t107 = t78 * pkin(6);
t83 = sin(qJ(4));
t106 = t78 * t83;
t105 = t78 * t84;
t86 = cos(qJ(4));
t104 = t78 * t86;
t81 = cos(pkin(5));
t103 = t81 * t84;
t102 = t81 * t87;
t64 = pkin(2) * t103 + (-pkin(6) - qJ(3)) * t78;
t73 = pkin(2) * t87 + pkin(1);
t77 = sin(pkin(9));
t80 = cos(pkin(9));
t99 = t80 * t64 + t77 * t73;
t98 = t81 * pkin(6) + qJ(1);
t97 = -t64 * t77 + t80 * t73;
t96 = pkin(2) * t105 + t81 * qJ(3) + t98;
t95 = g(1) * t77 - g(2) * t80;
t94 = t76 * t87 + t79 * t84;
t92 = t66 * t81;
t51 = -t77 * t94 + t80 * t92;
t63 = t94 * t81;
t52 = t63 * t80 + t66 * t77;
t93 = t52 * pkin(3) - pkin(7) * t51 + t99;
t53 = -t77 * t92 - t80 * t94;
t54 = -t63 * t77 + t66 * t80;
t91 = t54 * pkin(3) - pkin(7) * t53 + t97;
t45 = t104 * t80 + t52 * t83;
t47 = -t104 * t77 + t54 * t83;
t62 = t94 * t78;
t55 = t62 * t83 - t81 * t86;
t90 = g(1) * t47 + g(2) * t45 + g(3) * t55;
t61 = t66 * t78;
t89 = g(1) * t53 + g(2) * t51 + g(3) * t61;
t88 = t62 * pkin(3) - pkin(7) * t61 + t96;
t85 = cos(qJ(5));
t82 = sin(qJ(5));
t60 = -g(3) * t81 - t78 * t95;
t56 = t62 * t86 + t81 * t83;
t48 = t106 * t77 + t54 * t86;
t46 = -t106 * t80 + t52 * t86;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t77, t95, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t103 * t77 + t80 * t87) - g(2) * (t103 * t80 + t77 * t87) - g(3) * t105, -g(1) * (-t102 * t77 - t80 * t84) - g(2) * (t102 * t80 - t77 * t84) - g(3) * t78 * t87, t60, -g(1) * (pkin(1) * t80 + t107 * t77) - g(2) * (pkin(1) * t77 - t107 * t80) - g(3) * t98, 0, 0, 0, 0, 0, 0, -g(1) * t54 - g(2) * t52 - g(3) * t62, -t89, t60, -g(1) * t97 - g(2) * t99 - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * t48 - g(2) * t46 - g(3) * t56, t90, t89, -g(1) * t91 - g(2) * t93 - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (t48 * t85 - t53 * t82) - g(2) * (t46 * t85 - t51 * t82) - g(3) * (t56 * t85 - t61 * t82), -g(1) * (-t48 * t82 - t53 * t85) - g(2) * (-t46 * t82 - t51 * t85) - g(3) * (-t56 * t82 - t61 * t85), -t90, -g(1) * (pkin(4) * t48 + pkin(8) * t47 + t91) - g(2) * (pkin(4) * t46 + pkin(8) * t45 + t93) - g(3) * (pkin(4) * t56 + pkin(8) * t55 + t88);];
U_reg = t1;

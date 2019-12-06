% Calculate inertial parameters regressor of potential energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:21
% EndTime: 2019-12-05 16:23:21
% DurationCPUTime: 0.17s
% Computational Cost: add. (132->75), mult. (169->100), div. (0->0), fcn. (168->10), ass. (0->35)
t68 = sin(qJ(2));
t84 = g(3) * t68;
t67 = sin(qJ(3));
t83 = t67 * pkin(3);
t69 = cos(qJ(3));
t53 = t69 * pkin(3) + pkin(2);
t82 = g(3) * qJ(1);
t66 = -qJ(4) - pkin(6);
t62 = -pkin(7) + t66;
t81 = t62 * t68;
t64 = sin(pkin(8));
t80 = t64 * t67;
t70 = cos(qJ(2));
t79 = t64 * t70;
t65 = cos(pkin(8));
t78 = t65 * t70;
t77 = t66 * t68;
t76 = t67 * t70;
t75 = t69 * t70;
t74 = t65 * pkin(1) + t64 * pkin(5);
t63 = qJ(3) + pkin(9);
t58 = t64 * pkin(1);
t73 = -t65 * pkin(5) + t58;
t72 = pkin(2) * t70 + pkin(6) * t68;
t71 = g(1) * t65 + g(2) * t64;
t56 = qJ(5) + t63;
t55 = cos(t63);
t54 = sin(t63);
t52 = cos(t56);
t51 = sin(t56);
t50 = g(1) * t64 - g(2) * t65;
t49 = pkin(4) * t54 + t83;
t48 = pkin(4) * t55 + t53;
t47 = -g(3) * t70 + t71 * t68;
t1 = [0, 0, 0, 0, 0, 0, -t71, t50, -g(3), -t82, 0, 0, 0, 0, 0, 0, -t71 * t70 - t84, t47, -t50, -g(1) * t74 - g(2) * t73 - t82, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t75 + t80) - g(2) * (t64 * t75 - t65 * t67) - t69 * t84, -g(1) * (t64 * t69 - t65 * t76) - g(2) * (-t64 * t76 - t65 * t69) + t67 * t84, -t47, -g(1) * (t72 * t65 + t74) - g(2) * (t72 * t64 + t73) - g(3) * (t68 * pkin(2) - t70 * pkin(6) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t54 + t55 * t78) - g(2) * (-t65 * t54 + t55 * t79) - t55 * t84, -g(1) * (-t54 * t78 + t64 * t55) - g(2) * (-t54 * t79 - t65 * t55) + t54 * t84, -t47, -g(1) * (pkin(3) * t80 + t74) - g(2) * (t53 * t79 - t64 * t77 + t58) - g(3) * (t68 * t53 + t70 * t66 + qJ(1)) + (-g(1) * (t53 * t70 - t77) - g(2) * (-pkin(5) - t83)) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t51 + t52 * t78) - g(2) * (-t65 * t51 + t52 * t79) - t52 * t84, -g(1) * (-t51 * t78 + t64 * t52) - g(2) * (-t51 * t79 - t65 * t52) + t51 * t84, -t47, -g(1) * (t64 * t49 + t74) - g(2) * (t48 * t79 - t64 * t81 + t58) - g(3) * (t68 * t48 + t70 * t62 + qJ(1)) + (-g(1) * (t48 * t70 - t81) - g(2) * (-pkin(5) - t49)) * t65;];
U_reg = t1;

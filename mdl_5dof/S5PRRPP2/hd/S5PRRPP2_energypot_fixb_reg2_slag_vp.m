% Calculate inertial parameters regressor of potential energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:08
% EndTime: 2019-12-05 16:10:08
% DurationCPUTime: 0.12s
% Computational Cost: add. (129->61), mult. (184->80), div. (0->0), fcn. (187->8), ass. (0->36)
t67 = sin(qJ(2));
t85 = g(3) * t67;
t84 = g(3) * qJ(1);
t63 = sin(pkin(7));
t66 = sin(qJ(3));
t83 = t63 * t66;
t69 = cos(qJ(2));
t82 = t63 * t69;
t64 = cos(pkin(7));
t81 = t64 * t69;
t65 = -qJ(4) - pkin(6);
t80 = t65 * t67;
t79 = t66 * t69;
t68 = cos(qJ(3));
t78 = t68 * t69;
t77 = t64 * pkin(1) + t63 * pkin(5);
t55 = t68 * pkin(3) + pkin(2);
t76 = t67 * t55 + t69 * t65 + qJ(1);
t59 = t63 * pkin(1);
t75 = -t64 * pkin(5) + t59;
t74 = pkin(2) * t69 + pkin(6) * t67;
t73 = g(1) * t64 + g(2) * t63;
t72 = pkin(3) * t83 + t55 * t81 - t64 * t80 + t77;
t62 = qJ(3) + pkin(8);
t56 = sin(t62);
t57 = cos(t62);
t43 = t56 * t82 + t64 * t57;
t45 = t56 * t81 - t63 * t57;
t71 = g(1) * t45 + g(2) * t43 + t56 * t85;
t70 = -t63 * t80 + t55 * t82 + t59 + (-pkin(3) * t66 - pkin(5)) * t64;
t50 = g(1) * t63 - g(2) * t64;
t47 = -g(3) * t69 + t73 * t67;
t46 = t63 * t56 + t57 * t81;
t44 = -t64 * t56 + t57 * t82;
t42 = -g(1) * t46 - g(2) * t44 - t57 * t85;
t1 = [0, 0, 0, 0, 0, 0, -t73, t50, -g(3), -t84, 0, 0, 0, 0, 0, 0, -t73 * t69 - t85, t47, -t50, -g(1) * t77 - g(2) * t75 - t84, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t78 + t83) - g(2) * (t63 * t78 - t64 * t66) - t68 * t85, -g(1) * (t63 * t68 - t64 * t79) - g(2) * (-t63 * t79 - t64 * t68) + t66 * t85, -t47, -g(1) * (t74 * t64 + t77) - g(2) * (t74 * t63 + t75) - g(3) * (t67 * pkin(2) - t69 * pkin(6) + qJ(1)), 0, 0, 0, 0, 0, 0, t42, t71, -t47, -g(1) * t72 - g(2) * t70 - g(3) * t76, 0, 0, 0, 0, 0, 0, t42, -t47, -t71, -g(1) * (t46 * pkin(4) + t45 * qJ(5) + t72) - g(2) * (t44 * pkin(4) + t43 * qJ(5) + t70) - g(3) * ((pkin(4) * t57 + qJ(5) * t56) * t67 + t76);];
U_reg = t1;

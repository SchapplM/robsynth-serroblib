% Calculate inertial parameters regressor of potential energy for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:19
% EndTime: 2019-12-05 16:46:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (129->49), mult. (156->61), div. (0->0), fcn. (155->8), ass. (0->34)
t60 = qJ(2) + qJ(3);
t56 = sin(t60);
t57 = cos(t60);
t84 = pkin(3) * t57 + pkin(7) * t56;
t81 = g(3) * t56;
t80 = g(3) * qJ(1);
t61 = sin(pkin(8));
t63 = sin(qJ(4));
t79 = t61 * t63;
t65 = cos(qJ(4));
t78 = t61 * t65;
t62 = cos(pkin(8));
t77 = t62 * t63;
t76 = t62 * t65;
t66 = cos(qJ(2));
t55 = t66 * pkin(2) + pkin(1);
t67 = -pkin(6) - pkin(5);
t75 = t61 * t55 + t62 * t67;
t64 = sin(qJ(2));
t74 = t64 * pkin(2) + qJ(1);
t73 = t62 * t55 - t61 * t67;
t72 = t84 * t61 + t75;
t71 = g(1) * t62 + g(2) * t61;
t70 = t56 * pkin(3) - t57 * pkin(7) + t74;
t69 = t84 * t62 + t73;
t40 = t57 * t79 + t76;
t42 = t57 * t77 - t78;
t68 = g(1) * t42 + g(2) * t40 + t63 * t81;
t44 = g(1) * t61 - g(2) * t62;
t43 = t57 * t76 + t79;
t41 = t57 * t78 - t77;
t39 = -g(3) * t57 + t71 * t56;
t38 = -g(1) * t43 - g(2) * t41 - t65 * t81;
t1 = [0, 0, 0, 0, 0, 0, -t71, t44, -g(3), -t80, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t71 * t66, -g(3) * t66 + t71 * t64, -t44, -g(1) * (t62 * pkin(1) + t61 * pkin(5)) - g(2) * (t61 * pkin(1) - t62 * pkin(5)) - t80, 0, 0, 0, 0, 0, 0, -t71 * t57 - t81, t39, -t44, -g(1) * t73 - g(2) * t75 - g(3) * t74, 0, 0, 0, 0, 0, 0, t38, t68, -t39, -g(1) * t69 - g(2) * t72 - g(3) * t70, 0, 0, 0, 0, 0, 0, t38, -t39, -t68, -g(1) * (t43 * pkin(4) + t42 * qJ(5) + t69) - g(2) * (t41 * pkin(4) + t40 * qJ(5) + t72) - g(3) * ((pkin(4) * t65 + qJ(5) * t63) * t56 + t70);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5PRRRP6
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:16
% EndTime: 2019-12-05 16:52:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (129->61), mult. (184->80), div. (0->0), fcn. (187->8), ass. (0->36)
t67 = sin(qJ(2));
t86 = g(3) * t67;
t85 = g(3) * qJ(1);
t64 = sin(pkin(8));
t66 = sin(qJ(3));
t84 = t64 * t66;
t69 = cos(qJ(2));
t83 = t64 * t69;
t65 = cos(pkin(8));
t82 = t65 * t69;
t81 = t66 * t69;
t70 = -pkin(7) - pkin(6);
t80 = t67 * t70;
t68 = cos(qJ(3));
t79 = t68 * t69;
t78 = t65 * pkin(1) + t64 * pkin(5);
t55 = t68 * pkin(3) + pkin(2);
t77 = t67 * t55 + t69 * t70 + qJ(1);
t60 = t64 * pkin(1);
t76 = -t65 * pkin(5) + t60;
t75 = pkin(2) * t69 + pkin(6) * t67;
t74 = g(1) * t65 + g(2) * t64;
t73 = pkin(3) * t84 + t55 * t82 - t65 * t80 + t78;
t63 = qJ(3) + qJ(4);
t57 = sin(t63);
t58 = cos(t63);
t44 = t57 * t83 + t65 * t58;
t46 = t57 * t82 - t64 * t58;
t72 = g(1) * t46 + g(2) * t44 + t57 * t86;
t71 = -t64 * t80 + t55 * t83 + t60 + (-pkin(3) * t66 - pkin(5)) * t65;
t51 = g(1) * t64 - g(2) * t65;
t48 = -g(3) * t69 + t74 * t67;
t47 = t64 * t57 + t58 * t82;
t45 = -t65 * t57 + t58 * t83;
t43 = -g(1) * t47 - g(2) * t45 - t58 * t86;
t1 = [0, 0, 0, 0, 0, 0, -t74, t51, -g(3), -t85, 0, 0, 0, 0, 0, 0, -t74 * t69 - t86, t48, -t51, -g(1) * t78 - g(2) * t76 - t85, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t79 + t84) - g(2) * (t64 * t79 - t65 * t66) - t68 * t86, -g(1) * (t64 * t68 - t65 * t81) - g(2) * (-t64 * t81 - t65 * t68) + t66 * t86, -t48, -g(1) * (t75 * t65 + t78) - g(2) * (t75 * t64 + t76) - g(3) * (t67 * pkin(2) - t69 * pkin(6) + qJ(1)), 0, 0, 0, 0, 0, 0, t43, t72, -t48, -g(1) * t73 - g(2) * t71 - g(3) * t77, 0, 0, 0, 0, 0, 0, t43, -t48, -t72, -g(1) * (t47 * pkin(4) + t46 * qJ(5) + t73) - g(2) * (t45 * pkin(4) + t44 * qJ(5) + t71) - g(3) * ((pkin(4) * t58 + qJ(5) * t57) * t67 + t77);];
U_reg = t1;

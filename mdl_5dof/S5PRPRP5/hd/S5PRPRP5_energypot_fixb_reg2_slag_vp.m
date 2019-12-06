% Calculate inertial parameters regressor of potential energy for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:40
% EndTime: 2019-12-05 15:38:40
% DurationCPUTime: 0.13s
% Computational Cost: add. (129->61), mult. (184->78), div. (0->0), fcn. (187->8), ass. (0->34)
t66 = sin(qJ(2));
t81 = g(3) * t66;
t80 = g(3) * qJ(1);
t61 = sin(pkin(8));
t62 = sin(pkin(7));
t79 = t62 * t61;
t67 = cos(qJ(2));
t78 = t62 * t67;
t64 = cos(pkin(7));
t77 = t64 * t67;
t65 = -pkin(6) - qJ(3);
t76 = t65 * t66;
t75 = t64 * pkin(1) + t62 * pkin(5);
t63 = cos(pkin(8));
t52 = t63 * pkin(3) + pkin(2);
t74 = t66 * t52 + t67 * t65 + qJ(1);
t57 = t62 * pkin(1);
t73 = -t64 * pkin(5) + t57;
t72 = g(1) * t64 + g(2) * t62;
t71 = pkin(2) * t67 + qJ(3) * t66;
t70 = pkin(3) * t79 + t52 * t77 - t64 * t76 + t75;
t60 = pkin(8) + qJ(4);
t54 = sin(t60);
t55 = cos(t60);
t41 = t54 * t78 + t64 * t55;
t43 = t54 * t77 - t62 * t55;
t69 = g(1) * t43 + g(2) * t41 + t54 * t81;
t68 = -t62 * t76 + t52 * t78 + t57 + (-pkin(3) * t61 - pkin(5)) * t64;
t48 = g(1) * t62 - g(2) * t64;
t45 = -g(3) * t67 + t66 * t72;
t44 = t62 * t54 + t55 * t77;
t42 = -t64 * t54 + t55 * t78;
t40 = -g(1) * t44 - g(2) * t42 - t55 * t81;
t1 = [0, 0, 0, 0, 0, 0, -t72, t48, -g(3), -t80, 0, 0, 0, 0, 0, 0, -t67 * t72 - t81, t45, -t48, -g(1) * t75 - g(2) * t73 - t80, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t77 + t79) - g(2) * (-t64 * t61 + t63 * t78) - t63 * t81, -g(1) * (-t61 * t77 + t62 * t63) - g(2) * (-t61 * t78 - t64 * t63) + t61 * t81, -t45, -g(1) * (t64 * t71 + t75) - g(2) * (t62 * t71 + t73) - g(3) * (t66 * pkin(2) - t67 * qJ(3) + qJ(1)), 0, 0, 0, 0, 0, 0, t40, t69, -t45, -g(1) * t70 - g(2) * t68 - g(3) * t74, 0, 0, 0, 0, 0, 0, t40, -t45, -t69, -g(1) * (t44 * pkin(4) + t43 * qJ(5) + t70) - g(2) * (t42 * pkin(4) + t41 * qJ(5) + t68) - g(3) * ((pkin(4) * t55 + qJ(5) * t54) * t66 + t74);];
U_reg = t1;

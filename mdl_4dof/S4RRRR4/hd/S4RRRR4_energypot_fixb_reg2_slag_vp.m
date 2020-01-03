% Calculate inertial parameters regressor of potential energy for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:14
% EndTime: 2019-12-31 17:26:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (73->39), mult. (91->51), div. (0->0), fcn. (85->8), ass. (0->25)
t69 = g(3) * pkin(4);
t51 = qJ(2) + qJ(3);
t47 = sin(t51);
t68 = g(3) * t47;
t53 = sin(qJ(2));
t67 = t53 * pkin(2) + pkin(4);
t52 = sin(qJ(4));
t54 = sin(qJ(1));
t66 = t54 * t52;
t55 = cos(qJ(4));
t65 = t54 * t55;
t57 = cos(qJ(1));
t64 = t57 * t52;
t63 = t57 * t55;
t56 = cos(qJ(2));
t45 = t56 * pkin(2) + pkin(1);
t58 = -pkin(6) - pkin(5);
t62 = t54 * t45 + t57 * t58;
t61 = t57 * t45 - t54 * t58;
t48 = cos(t51);
t60 = pkin(3) * t48 + pkin(7) * t47;
t59 = g(1) * t57 + g(2) * t54;
t42 = g(1) * t54 - g(2) * t57;
t41 = -g(3) * t48 + t59 * t47;
t1 = [0, 0, 0, 0, 0, 0, -t59, t42, -g(3), -t69, 0, 0, 0, 0, 0, 0, -g(3) * t53 - t59 * t56, -g(3) * t56 + t59 * t53, -t42, -g(1) * (t57 * pkin(1) + t54 * pkin(5)) - g(2) * (t54 * pkin(1) - t57 * pkin(5)) - t69, 0, 0, 0, 0, 0, 0, -t59 * t48 - t68, t41, -t42, -g(1) * t61 - g(2) * t62 - g(3) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (t48 * t63 + t66) - g(2) * (t48 * t65 - t64) - t55 * t68, -g(1) * (-t48 * t64 + t65) - g(2) * (-t48 * t66 - t63) + t52 * t68, -t41, -g(1) * (t60 * t57 + t61) - g(2) * (t60 * t54 + t62) - g(3) * (t47 * pkin(3) - t48 * pkin(7) + t67);];
U_reg = t1;

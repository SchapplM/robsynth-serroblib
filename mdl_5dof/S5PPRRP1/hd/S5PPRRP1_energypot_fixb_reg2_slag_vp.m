% Calculate inertial parameters regressor of potential energy for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:13
% EndTime: 2019-12-05 15:07:13
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->52), mult. (143->65), div. (0->0), fcn. (138->8), ass. (0->31)
t60 = cos(qJ(4));
t47 = t60 * pkin(4) + pkin(3);
t52 = pkin(8) + qJ(3);
t48 = sin(t52);
t49 = cos(t52);
t57 = -qJ(5) - pkin(6);
t74 = t47 * t49 - t48 * t57;
t73 = g(3) * t48;
t72 = g(3) * qJ(1);
t54 = sin(pkin(7));
t59 = sin(qJ(4));
t69 = t54 * t59;
t68 = t54 * t60;
t56 = cos(pkin(7));
t67 = t56 * t59;
t66 = t56 * t60;
t55 = cos(pkin(8));
t46 = t55 * pkin(2) + pkin(1);
t58 = -pkin(5) - qJ(2);
t65 = t54 * t46 + t56 * t58;
t53 = sin(pkin(8));
t64 = t53 * pkin(2) + qJ(1);
t43 = t56 * t46;
t63 = -t54 * t58 + t43;
t62 = pkin(3) * t49 + pkin(6) * t48;
t61 = g(1) * t56 + g(2) * t54;
t41 = g(1) * t54 - g(2) * t56;
t40 = -g(3) * t49 + t61 * t48;
t39 = -g(1) * (t49 * t66 + t69) - g(2) * (t49 * t68 - t67) - t60 * t73;
t38 = -g(1) * (-t49 * t67 + t68) - g(2) * (-t49 * t69 - t66) + t59 * t73;
t1 = [0, 0, 0, 0, 0, 0, -t61, t41, -g(3), -t72, 0, 0, 0, 0, 0, 0, -g(3) * t53 - t61 * t55, -g(3) * t55 + t61 * t53, -t41, -g(1) * (t56 * pkin(1) + t54 * qJ(2)) - g(2) * (t54 * pkin(1) - t56 * qJ(2)) - t72, 0, 0, 0, 0, 0, 0, -t61 * t49 - t73, t40, -t41, -g(1) * t63 - g(2) * t65 - g(3) * t64, 0, 0, 0, 0, 0, 0, t39, t38, -t40, -g(1) * (t62 * t56 + t63) - g(2) * (t62 * t54 + t65) - g(3) * (t48 * pkin(3) - t49 * pkin(6) + t64), 0, 0, 0, 0, 0, 0, t39, t38, -t40, -g(1) * (t74 * t56 + t43) - g(2) * (-pkin(4) * t67 + t65) - g(3) * (t48 * t47 + t49 * t57 + t64) + (-g(1) * (pkin(4) * t59 - t58) - g(2) * t74) * t54;];
U_reg = t1;

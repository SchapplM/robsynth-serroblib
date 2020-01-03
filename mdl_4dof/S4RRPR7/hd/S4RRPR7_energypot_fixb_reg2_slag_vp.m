% Calculate inertial parameters regressor of potential energy for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:43
% EndTime: 2019-12-31 17:06:43
% DurationCPUTime: 0.06s
% Computational Cost: add. (73->39), mult. (91->51), div. (0->0), fcn. (85->8), ass. (0->25)
t69 = g(3) * pkin(4);
t51 = qJ(2) + pkin(7);
t47 = sin(t51);
t68 = g(3) * t47;
t54 = sin(qJ(2));
t67 = t54 * pkin(2) + pkin(4);
t53 = sin(qJ(4));
t55 = sin(qJ(1));
t66 = t55 * t53;
t56 = cos(qJ(4));
t65 = t55 * t56;
t58 = cos(qJ(1));
t64 = t58 * t53;
t63 = t58 * t56;
t57 = cos(qJ(2));
t46 = t57 * pkin(2) + pkin(1);
t52 = -pkin(5) - qJ(3);
t62 = t55 * t46 + t58 * t52;
t61 = t58 * t46 - t55 * t52;
t48 = cos(t51);
t60 = pkin(3) * t48 + pkin(6) * t47;
t59 = g(1) * t58 + g(2) * t55;
t42 = g(1) * t55 - g(2) * t58;
t41 = -g(3) * t48 + t59 * t47;
t1 = [0, 0, 0, 0, 0, 0, -t59, t42, -g(3), -t69, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t59 * t57, -g(3) * t57 + t59 * t54, -t42, -g(1) * (t58 * pkin(1) + t55 * pkin(5)) - g(2) * (t55 * pkin(1) - t58 * pkin(5)) - t69, 0, 0, 0, 0, 0, 0, -t59 * t48 - t68, t41, -t42, -g(1) * t61 - g(2) * t62 - g(3) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (t48 * t63 + t66) - g(2) * (t48 * t65 - t64) - t56 * t68, -g(1) * (-t48 * t64 + t65) - g(2) * (-t48 * t66 - t63) + t53 * t68, -t41, -g(1) * (t60 * t58 + t61) - g(2) * (t60 * t55 + t62) - g(3) * (t47 * pkin(3) - t48 * pkin(6) + t67);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S4RRRR5
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:12
% EndTime: 2019-12-31 17:28:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (71->51), mult. (115->68), div. (0->0), fcn. (113->8), ass. (0->28)
t73 = g(3) * pkin(4);
t55 = sin(qJ(2));
t72 = g(3) * t55;
t60 = -pkin(7) - pkin(6);
t71 = t55 * t60;
t54 = sin(qJ(3));
t56 = sin(qJ(1));
t70 = t56 * t54;
t58 = cos(qJ(2));
t69 = t56 * t58;
t53 = qJ(3) + qJ(4);
t47 = sin(t53);
t59 = cos(qJ(1));
t68 = t59 * t47;
t48 = cos(t53);
t67 = t59 * t48;
t66 = t59 * t54;
t57 = cos(qJ(3));
t65 = t59 * t57;
t64 = t59 * pkin(1) + t56 * pkin(5);
t50 = t56 * pkin(1);
t63 = -t59 * pkin(5) + t50;
t62 = pkin(2) * t58 + pkin(6) * t55;
t61 = g(1) * t59 + g(2) * t56;
t46 = t57 * pkin(3) + pkin(2);
t45 = g(1) * t56 - g(2) * t59;
t44 = -g(3) * t58 + t61 * t55;
t1 = [0, 0, 0, 0, 0, 0, -t61, t45, -g(3), -t73, 0, 0, 0, 0, 0, 0, -t61 * t58 - t72, t44, -t45, -g(1) * t64 - g(2) * t63 - t73, 0, 0, 0, 0, 0, 0, -g(1) * (t58 * t65 + t70) - g(2) * (t57 * t69 - t66) - t57 * t72, -g(1) * (t56 * t57 - t58 * t66) - g(2) * (-t54 * t69 - t65) + t54 * t72, -t44, -g(1) * (t62 * t59 + t64) - g(2) * (t62 * t56 + t63) - g(3) * (t55 * pkin(2) - t58 * pkin(6) + pkin(4)), 0, 0, 0, 0, 0, 0, -g(1) * (t56 * t47 + t58 * t67) - g(2) * (t48 * t69 - t68) - t48 * t72, -g(1) * (t56 * t48 - t58 * t68) - g(2) * (-t47 * t69 - t67) + t47 * t72, -t44, -g(1) * (pkin(3) * t70 + t64) - g(2) * (t46 * t69 - t56 * t71 + t50) - g(3) * (t55 * t46 + t58 * t60 + pkin(4)) + (-g(1) * (t46 * t58 - t71) - g(2) * (-pkin(3) * t54 - pkin(5))) * t59;];
U_reg = t1;

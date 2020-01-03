% Calculate inertial parameters regressor of potential energy for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:00
% EndTime: 2019-12-31 17:44:00
% DurationCPUTime: 0.11s
% Computational Cost: add. (119->43), mult. (120->51), div. (0->0), fcn. (113->8), ass. (0->28)
t53 = qJ(1) + pkin(7);
t48 = sin(t53);
t54 = sin(pkin(8));
t69 = qJ(4) * t54;
t55 = cos(pkin(8));
t72 = t48 * t55;
t76 = pkin(3) * t72 + t48 * t69;
t49 = cos(t53);
t65 = g(1) * t49 + g(2) * t48;
t56 = qJ(2) + pkin(5);
t73 = g(3) * t56;
t71 = t49 * t55;
t58 = sin(qJ(1));
t70 = t58 * pkin(1) + t48 * pkin(2);
t60 = cos(qJ(1));
t68 = t60 * pkin(1) + t49 * pkin(2) + t48 * qJ(3);
t67 = -t49 * qJ(3) + t70;
t66 = pkin(3) * t71 + t49 * t69 + t68;
t64 = -g(1) * t60 - g(2) * t58;
t57 = sin(qJ(5));
t59 = cos(qJ(5));
t63 = t54 * t59 - t55 * t57;
t62 = t54 * t57 + t55 * t59;
t61 = t54 * pkin(3) - t55 * qJ(4) + t56;
t39 = g(1) * t48 - g(2) * t49;
t38 = -g(3) * t54 - t55 * t65;
t37 = -g(3) * t55 + t54 * t65;
t1 = [0, 0, 0, 0, 0, 0, t64, g(1) * t58 - g(2) * t60, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t65, t39, -g(3), pkin(1) * t64 - t73, 0, 0, 0, 0, 0, 0, t38, t37, -t39, -g(1) * t68 - g(2) * t67 - t73, 0, 0, 0, 0, 0, 0, t38, -t39, -t37, -g(1) * t66 - g(2) * (t67 + t76) - g(3) * t61, 0, 0, 0, 0, 0, 0, -g(3) * t63 - t65 * t62, g(3) * t62 - t65 * t63, t39, -g(1) * (pkin(4) * t71 - t48 * pkin(6) + t66) - g(2) * (pkin(4) * t72 + (pkin(6) - qJ(3)) * t49 + t70 + t76) - g(3) * (t54 * pkin(4) + t61);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:54
% EndTime: 2019-12-31 17:36:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (119->43), mult. (120->51), div. (0->0), fcn. (113->8), ass. (0->28)
t49 = pkin(7) + qJ(2);
t44 = sin(t49);
t50 = sin(pkin(8));
t65 = qJ(4) * t50;
t52 = cos(pkin(8));
t68 = t44 * t52;
t72 = pkin(3) * t68 + t44 * t65;
t45 = cos(t49);
t61 = g(1) * t45 + g(2) * t44;
t54 = pkin(5) + qJ(1);
t69 = g(3) * t54;
t67 = t45 * t52;
t51 = sin(pkin(7));
t66 = t51 * pkin(1) + t44 * pkin(2);
t53 = cos(pkin(7));
t64 = t53 * pkin(1) + t45 * pkin(2) + t44 * qJ(3);
t63 = -t45 * qJ(3) + t66;
t62 = pkin(3) * t67 + t45 * t65 + t64;
t60 = -g(1) * t53 - g(2) * t51;
t55 = sin(qJ(5));
t56 = cos(qJ(5));
t59 = t50 * t56 - t52 * t55;
t58 = t50 * t55 + t52 * t56;
t57 = t50 * pkin(3) - t52 * qJ(4) + t54;
t35 = g(1) * t44 - g(2) * t45;
t34 = -g(3) * t50 - t52 * t61;
t33 = -g(3) * t52 + t50 * t61;
t1 = [0, 0, 0, 0, 0, 0, t60, g(1) * t51 - g(2) * t53, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t61, t35, -g(3), pkin(1) * t60 - t69, 0, 0, 0, 0, 0, 0, t34, t33, -t35, -g(1) * t64 - g(2) * t63 - t69, 0, 0, 0, 0, 0, 0, t34, -t35, -t33, -g(1) * t62 - g(2) * (t63 + t72) - g(3) * t57, 0, 0, 0, 0, 0, 0, -g(3) * t59 - t61 * t58, g(3) * t58 - t61 * t59, t35, -g(1) * (pkin(4) * t67 - t44 * pkin(6) + t62) - g(2) * (pkin(4) * t68 + (pkin(6) - qJ(3)) * t45 + t66 + t72) - g(3) * (t50 * pkin(4) + t57);];
U_reg = t1;

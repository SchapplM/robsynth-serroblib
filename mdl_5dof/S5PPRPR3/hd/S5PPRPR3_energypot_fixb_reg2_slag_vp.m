% Calculate inertial parameters regressor of potential energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:21
% EndTime: 2019-12-05 15:05:21
% DurationCPUTime: 0.19s
% Computational Cost: add. (142->70), mult. (205->97), div. (0->0), fcn. (216->10), ass. (0->41)
t50 = sin(pkin(8));
t78 = g(3) * t50;
t77 = g(3) * qJ(1);
t54 = -qJ(4) - pkin(5);
t76 = t50 * t54;
t55 = sin(qJ(5));
t75 = t50 * t55;
t57 = cos(qJ(5));
t74 = t50 * t57;
t51 = sin(pkin(7));
t52 = cos(pkin(8));
t73 = t51 * t52;
t56 = sin(qJ(3));
t72 = t51 * t56;
t58 = cos(qJ(3));
t71 = t51 * t58;
t49 = qJ(3) + pkin(9);
t43 = sin(t49);
t53 = cos(pkin(7));
t70 = t53 * t43;
t44 = cos(t49);
t69 = t53 * t44;
t68 = t53 * t56;
t67 = t53 * t58;
t66 = t53 * pkin(1) + t51 * qJ(2);
t42 = t58 * pkin(3) + pkin(2);
t65 = t50 * t42 + t52 * t54 + qJ(1);
t46 = t51 * pkin(1);
t64 = -t53 * qJ(2) + t46;
t63 = pkin(2) * t52 + pkin(5) * t50;
t62 = g(1) * t53 + g(2) * t51;
t61 = pkin(3) * t72 + t66 + (t42 * t52 - t76) * t53;
t30 = t43 * t73 + t69;
t32 = -t51 * t44 + t52 * t70;
t60 = g(1) * t32 + g(2) * t30 + t43 * t78;
t59 = -t51 * t76 + t42 * t73 + t46 + (-pkin(3) * t56 - qJ(2)) * t53;
t37 = g(1) * t51 - g(2) * t53;
t34 = -g(3) * t52 + t62 * t50;
t33 = t51 * t43 + t52 * t69;
t31 = t44 * t73 - t70;
t1 = [0, 0, 0, 0, 0, 0, -t62, t37, -g(3), -t77, 0, 0, 0, 0, 0, 0, -t62 * t52 - t78, t34, -t37, -g(1) * t66 - g(2) * t64 - t77, 0, 0, 0, 0, 0, 0, -g(1) * (t52 * t67 + t72) - g(2) * (t52 * t71 - t68) - t58 * t78, -g(1) * (-t52 * t68 + t71) - g(2) * (-t52 * t72 - t67) + t56 * t78, -t34, -g(1) * (t63 * t53 + t66) - g(2) * (t63 * t51 + t64) - g(3) * (t50 * pkin(2) - t52 * pkin(5) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * t33 - g(2) * t31 - t44 * t78, t60, -t34, -g(1) * t61 - g(2) * t59 - g(3) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (t33 * t57 + t53 * t75) - g(2) * (t31 * t57 + t51 * t75) - g(3) * (t44 * t74 - t52 * t55), -g(1) * (-t33 * t55 + t53 * t74) - g(2) * (-t31 * t55 + t51 * t74) - g(3) * (-t44 * t75 - t52 * t57), -t60, -g(1) * (t33 * pkin(4) + t32 * pkin(6) + t61) - g(2) * (t31 * pkin(4) + t30 * pkin(6) + t59) - g(3) * ((pkin(4) * t44 + pkin(6) * t43) * t50 + t65);];
U_reg = t1;

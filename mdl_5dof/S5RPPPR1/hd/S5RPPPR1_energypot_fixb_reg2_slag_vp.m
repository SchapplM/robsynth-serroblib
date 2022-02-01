% Calculate inertial parameters regressor of potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:37
% EndTime: 2022-01-20 09:12:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (138->59), mult. (130->79), div. (0->0), fcn. (125->10), ass. (0->32)
t61 = sin(pkin(8));
t81 = g(3) * t61;
t64 = qJ(2) + pkin(5);
t80 = g(3) * t64;
t59 = qJ(1) + pkin(7);
t52 = sin(t59);
t60 = sin(pkin(9));
t79 = t52 * t60;
t63 = cos(pkin(8));
t78 = t52 * t63;
t54 = cos(t59);
t77 = t54 * t63;
t76 = t60 * t63;
t65 = -pkin(6) - qJ(4);
t75 = t61 * t65;
t62 = cos(pkin(9));
t74 = t62 * t63;
t66 = sin(qJ(1));
t73 = t66 * pkin(1) + t52 * pkin(2);
t67 = cos(qJ(1));
t72 = t67 * pkin(1) + t54 * pkin(2) + t52 * qJ(3);
t71 = -t54 * qJ(3) + t73;
t70 = g(1) * t54 + g(2) * t52;
t69 = -g(1) * t67 - g(2) * t66;
t68 = pkin(3) * t63 + qJ(4) * t61;
t58 = pkin(9) + qJ(5);
t53 = cos(t58);
t51 = sin(t58);
t50 = t62 * pkin(4) + pkin(3);
t46 = g(1) * t52 - g(2) * t54;
t45 = -g(3) * t63 + t70 * t61;
t1 = [0, 0, 0, 0, 0, 0, t69, g(1) * t66 - g(2) * t67, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t70, t46, -g(3), t69 * pkin(1) - t80, 0, 0, 0, 0, 0, 0, -t70 * t63 - t81, t45, -t46, -g(1) * t72 - g(2) * t71 - t80, 0, 0, 0, 0, 0, 0, -g(1) * (t54 * t74 + t79) - g(2) * (t52 * t74 - t54 * t60) - t62 * t81, -g(1) * (t52 * t62 - t54 * t76) - g(2) * (-t52 * t76 - t54 * t62) + t60 * t81, -t45, -g(1) * (t68 * t54 + t72) - g(2) * (t68 * t52 + t71) - g(3) * (t61 * pkin(3) - t63 * qJ(4) + t64), 0, 0, 0, 0, 0, 0, -g(1) * (t52 * t51 + t53 * t77) - g(2) * (-t54 * t51 + t53 * t78) - t53 * t81, -g(1) * (-t51 * t77 + t52 * t53) - g(2) * (-t51 * t78 - t54 * t53) + t51 * t81, -t45, -g(1) * (pkin(4) * t79 + t72) - g(2) * (t50 * t78 - t52 * t75 + t73) - g(3) * (t61 * t50 + t63 * t65 + t64) + (-g(1) * (t50 * t63 - t75) - g(2) * (-pkin(4) * t60 - qJ(3))) * t54;];
U_reg = t1;

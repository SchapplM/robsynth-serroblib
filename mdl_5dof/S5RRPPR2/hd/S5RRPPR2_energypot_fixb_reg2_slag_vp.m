% Calculate inertial parameters regressor of potential energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:35
% EndTime: 2020-01-03 11:57:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (134->46), mult. (94->59), div. (0->0), fcn. (85->10), ass. (0->28)
t55 = sin(pkin(9));
t56 = cos(pkin(9));
t74 = pkin(4) * t56 + pkin(7) * t55;
t73 = pkin(6) + pkin(5);
t53 = qJ(3) + t73;
t70 = g(1) * t53;
t69 = g(1) * t55;
t57 = sin(qJ(5));
t68 = t56 * t57;
t59 = cos(qJ(5));
t67 = t56 * t59;
t54 = qJ(1) + qJ(2);
t50 = sin(t54);
t58 = sin(qJ(1));
t66 = t58 * pkin(1) + pkin(2) * t50;
t49 = pkin(8) + t54;
t46 = sin(t49);
t65 = t46 * pkin(3) + t66;
t51 = cos(t54);
t60 = cos(qJ(1));
t64 = -t60 * pkin(1) - pkin(2) * t51;
t47 = cos(t49);
t63 = -g(2) * t46 + g(3) * t47;
t62 = -g(2) * t58 + g(3) * t60;
t61 = -t46 * qJ(4) + t64;
t44 = g(2) * t47 + g(3) * t46;
t43 = g(1) * t56 + t63 * t55;
t1 = [0, 0, 0, 0, 0, 0, t62, -g(2) * t60 - g(3) * t58, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -g(2) * t50 + g(3) * t51, -g(2) * t51 - g(3) * t50, -g(1), t62 * pkin(1) - g(1) * t73, 0, 0, 0, 0, 0, 0, t63, -t44, -g(1), -g(2) * t66 - g(3) * t64 - t70, 0, 0, 0, 0, 0, 0, t63 * t56 - t69, -t43, t44, -t70 - g(2) * (-t47 * qJ(4) + t65) - g(3) * (-t47 * pkin(3) + t61), 0, 0, 0, 0, 0, 0, -t59 * t69 - g(2) * (t46 * t67 - t47 * t57) - g(3) * (-t46 * t57 - t47 * t67), t57 * t69 - g(2) * (-t46 * t68 - t47 * t59) - g(3) * (-t46 * t59 + t47 * t68), t43, -g(1) * (t55 * pkin(4) - t56 * pkin(7) + t53) - g(2) * (t74 * t46 + t65) - g(3) * t61 + (g(2) * qJ(4) - g(3) * (-pkin(3) - t74)) * t47;];
U_reg = t1;

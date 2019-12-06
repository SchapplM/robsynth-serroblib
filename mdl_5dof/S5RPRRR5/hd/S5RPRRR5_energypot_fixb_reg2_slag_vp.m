% Calculate inertial parameters regressor of potential energy for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:27
% EndTime: 2019-12-05 18:16:27
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->41), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t67 = qJ(2) + pkin(5);
t55 = pkin(6) + t67;
t68 = g(1) * t55;
t56 = qJ(1) + pkin(9);
t50 = cos(t56);
t61 = cos(qJ(1));
t66 = t61 * pkin(1) + pkin(2) * t50;
t49 = sin(t56);
t59 = sin(qJ(1));
t65 = -t59 * pkin(1) - pkin(2) * t49;
t51 = qJ(3) + t56;
t46 = sin(t51);
t47 = cos(t51);
t64 = g(2) * t46 - g(3) * t47;
t63 = g(2) * t59 - g(3) * t61;
t62 = -pkin(8) - pkin(7);
t60 = cos(qJ(4));
t58 = sin(qJ(4));
t57 = qJ(4) + qJ(5);
t53 = cos(t57);
t52 = sin(t57);
t48 = t60 * pkin(4) + pkin(3);
t44 = g(2) * t47 + g(3) * t46;
t1 = [0, 0, 0, 0, 0, 0, t63, g(2) * t61 + g(3) * t59, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, g(2) * t49 - g(3) * t50, g(2) * t50 + g(3) * t49, -g(1), t63 * pkin(1) - g(1) * t67, 0, 0, 0, 0, 0, 0, t64, t44, -g(1), -g(2) * t65 - g(3) * t66 - t68, 0, 0, 0, 0, 0, 0, -g(1) * t58 + t64 * t60, -g(1) * t60 - t64 * t58, -t44, -t68 - g(2) * (-t46 * pkin(3) + t47 * pkin(7) + t65) - g(3) * (t47 * pkin(3) + t46 * pkin(7) + t66), 0, 0, 0, 0, 0, 0, -g(1) * t52 + t64 * t53, -g(1) * t53 - t64 * t52, -t44, -g(1) * (t58 * pkin(4) + t55) - g(2) * (-t46 * t48 - t47 * t62 + t65) - g(3) * (-t46 * t62 + t47 * t48 + t66);];
U_reg = t1;

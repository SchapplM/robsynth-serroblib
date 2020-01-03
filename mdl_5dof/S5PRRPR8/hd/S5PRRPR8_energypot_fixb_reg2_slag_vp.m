% Calculate inertial parameters regressor of potential energy for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:41
% EndTime: 2019-12-31 17:42:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (130->50), mult. (119->64), div. (0->0), fcn. (110->10), ass. (0->31)
t54 = -pkin(6) - pkin(5);
t47 = qJ(2) + qJ(3);
t40 = pkin(9) + t47;
t36 = sin(t40);
t66 = g(3) * t36;
t53 = cos(qJ(2));
t39 = t53 * pkin(2) + pkin(1);
t65 = g(3) * qJ(1);
t48 = sin(pkin(8));
t50 = sin(qJ(5));
t64 = t48 * t50;
t52 = cos(qJ(5));
t63 = t48 * t52;
t49 = cos(pkin(8));
t62 = t49 * t50;
t61 = t49 * t52;
t42 = cos(t47);
t33 = pkin(3) * t42 + t39;
t46 = -qJ(4) + t54;
t60 = t48 * t33 + t49 * t46;
t51 = sin(qJ(2));
t59 = t51 * pkin(2) + qJ(1);
t41 = sin(t47);
t58 = pkin(3) * t41 + t59;
t57 = t49 * t33 - t48 * t46;
t37 = cos(t40);
t56 = pkin(4) * t37 + pkin(7) * t36;
t55 = g(1) * t49 + g(2) * t48;
t34 = g(1) * t48 - g(2) * t49;
t30 = -g(3) * t37 + t55 * t36;
t1 = [0, 0, 0, 0, 0, 0, -t55, t34, -g(3), -t65, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t55 * t53, -g(3) * t53 + t55 * t51, -t34, -g(1) * (t49 * pkin(1) + t48 * pkin(5)) - g(2) * (t48 * pkin(1) - t49 * pkin(5)) - t65, 0, 0, 0, 0, 0, 0, -g(3) * t41 - t55 * t42, -g(3) * t42 + t55 * t41, -t34, -g(1) * (t49 * t39 - t48 * t54) - g(2) * (t48 * t39 + t49 * t54) - g(3) * t59, 0, 0, 0, 0, 0, 0, -t55 * t37 - t66, t30, -t34, -g(1) * t57 - g(2) * t60 - g(3) * t58, 0, 0, 0, 0, 0, 0, -g(1) * (t37 * t61 + t64) - g(2) * (t37 * t63 - t62) - t52 * t66, -g(1) * (-t37 * t62 + t63) - g(2) * (-t37 * t64 - t61) + t50 * t66, -t30, -g(1) * (t56 * t49 + t57) - g(2) * (t56 * t48 + t60) - g(3) * (t36 * pkin(4) - t37 * pkin(7) + t58);];
U_reg = t1;

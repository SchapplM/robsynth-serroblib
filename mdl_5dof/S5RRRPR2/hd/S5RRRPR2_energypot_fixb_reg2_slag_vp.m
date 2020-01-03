% Calculate inertial parameters regressor of potential energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:36
% EndTime: 2020-01-03 12:07:36
% DurationCPUTime: 0.05s
% Computational Cost: add. (117->38), mult. (63->41), div. (0->0), fcn. (50->10), ass. (0->24)
t59 = pkin(6) + pkin(5);
t56 = pkin(7) + t59;
t58 = g(1) * (qJ(4) + t56);
t46 = qJ(1) + qJ(2);
t41 = sin(t46);
t48 = sin(qJ(1));
t57 = t48 * pkin(1) + pkin(2) * t41;
t44 = qJ(3) + t46;
t39 = sin(t44);
t55 = pkin(3) * t39 + t57;
t42 = cos(t46);
t50 = cos(qJ(1));
t54 = -t50 * pkin(1) - pkin(2) * t42;
t38 = pkin(9) + t44;
t34 = sin(t38);
t35 = cos(t38);
t53 = g(2) * t34 - g(3) * t35;
t52 = -g(2) * t48 + g(3) * t50;
t40 = cos(t44);
t51 = -pkin(3) * t40 + t54;
t49 = cos(qJ(5));
t47 = sin(qJ(5));
t33 = g(2) * t35 + g(3) * t34;
t1 = [0, 0, 0, 0, 0, 0, t52, -g(2) * t50 - g(3) * t48, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -g(2) * t41 + g(3) * t42, -g(2) * t42 - g(3) * t41, -g(1), t52 * pkin(1) - g(1) * t59, 0, 0, 0, 0, 0, 0, -g(2) * t39 + g(3) * t40, -g(2) * t40 - g(3) * t39, -g(1), -g(1) * t56 - g(2) * t57 - g(3) * t54, 0, 0, 0, 0, 0, 0, -t53, -t33, -g(1), -g(2) * t55 - g(3) * t51 - t58, 0, 0, 0, 0, 0, 0, -g(1) * t47 - t53 * t49, -g(1) * t49 + t53 * t47, t33, -t58 - g(2) * (t34 * pkin(4) - t35 * pkin(8) + t55) - g(3) * (-t35 * pkin(4) - t34 * pkin(8) + t51);];
U_reg = t1;

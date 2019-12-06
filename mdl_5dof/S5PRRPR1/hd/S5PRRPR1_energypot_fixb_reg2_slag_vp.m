% Calculate inertial parameters regressor of potential energy for
% S5PRRPR1
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:49
% EndTime: 2019-12-05 16:15:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->42), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t61 = pkin(5) + qJ(1);
t49 = pkin(6) + t61;
t62 = g(3) * t49;
t51 = pkin(8) + qJ(2);
t43 = sin(t51);
t53 = sin(pkin(8));
t60 = t53 * pkin(1) + pkin(2) * t43;
t45 = cos(t51);
t55 = cos(pkin(8));
t59 = t55 * pkin(1) + pkin(2) * t45;
t46 = qJ(3) + t51;
t39 = sin(t46);
t40 = cos(t46);
t58 = g(1) * t40 + g(2) * t39;
t57 = -g(1) * t55 - g(2) * t53;
t56 = -pkin(7) - qJ(4);
t54 = cos(pkin(9));
t52 = sin(pkin(9));
t50 = pkin(9) + qJ(5);
t44 = cos(t50);
t42 = sin(t50);
t41 = t54 * pkin(4) + pkin(3);
t35 = g(1) * t39 - g(2) * t40;
t1 = [0, 0, 0, 0, 0, 0, t57, g(1) * t53 - g(2) * t55, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t43, g(1) * t43 - g(2) * t45, -g(3), t57 * pkin(1) - g(3) * t61, 0, 0, 0, 0, 0, 0, -t58, t35, -g(3), -g(1) * t59 - g(2) * t60 - t62, 0, 0, 0, 0, 0, 0, -g(3) * t52 - t58 * t54, -g(3) * t54 + t58 * t52, -t35, -g(1) * (t40 * pkin(3) + t39 * qJ(4) + t59) - g(2) * (t39 * pkin(3) - t40 * qJ(4) + t60) - t62, 0, 0, 0, 0, 0, 0, -g(3) * t42 - t58 * t44, -g(3) * t44 + t58 * t42, -t35, -g(1) * (-t39 * t56 + t40 * t41 + t59) - g(2) * (t39 * t41 + t40 * t56 + t60) - g(3) * (t52 * pkin(4) + t49);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:42
% EndTime: 2019-12-05 17:47:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (83->46), mult. (92->46), div. (0->0), fcn. (79->8), ass. (0->23)
t67 = g(3) * pkin(5);
t66 = pkin(2) + pkin(5);
t58 = sin(qJ(3));
t65 = t58 * pkin(3);
t57 = -qJ(4) - pkin(6);
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t64 = t61 * pkin(1) + t59 * qJ(2);
t56 = qJ(3) + pkin(8);
t60 = cos(qJ(3));
t63 = t60 * pkin(3) + t66;
t52 = t59 * pkin(1);
t62 = -t61 * qJ(2) + t52;
t44 = g(1) * t59 - g(2) * t61;
t55 = -pkin(7) + t57;
t50 = qJ(5) + t56;
t49 = cos(t56);
t48 = sin(t56);
t47 = cos(t50);
t46 = sin(t50);
t45 = g(1) * t61 + g(2) * t59;
t43 = pkin(4) * t48 + t65;
t1 = [0, 0, 0, 0, 0, 0, -t45, t44, -g(3), -t67, 0, 0, 0, 0, 0, 0, -g(3), t45, -t44, -g(1) * t64 - g(2) * t62 - t67, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t44 * t58, g(3) * t58 - t44 * t60, -t45, -g(1) * (t61 * pkin(6) + t64) - g(2) * (t59 * pkin(6) + t62) - g(3) * t66, 0, 0, 0, 0, 0, 0, -g(3) * t49 - t44 * t48, g(3) * t48 - t44 * t49, -t45, -g(1) * (-t61 * t57 + t59 * t65 + t64) - g(2) * (-t59 * t57 + t52 + (-qJ(2) - t65) * t61) - g(3) * t63, 0, 0, 0, 0, 0, 0, -g(3) * t47 - t44 * t46, g(3) * t46 - t44 * t47, -t45, -g(1) * (t59 * t43 - t61 * t55 + t64) - g(2) * (-t59 * t55 + t52 + (-qJ(2) - t43) * t61) - g(3) * (pkin(4) * t49 + t63);];
U_reg = t1;

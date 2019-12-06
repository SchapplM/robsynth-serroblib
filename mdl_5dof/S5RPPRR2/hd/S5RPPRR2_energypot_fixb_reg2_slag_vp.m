% Calculate inertial parameters regressor of potential energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:55
% EndTime: 2019-12-05 17:39:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (83->46), mult. (92->46), div. (0->0), fcn. (79->8), ass. (0->23)
t64 = g(3) * pkin(5);
t63 = pkin(2) + pkin(5);
t54 = sin(pkin(8));
t62 = t54 * pkin(3);
t56 = -pkin(6) - qJ(3);
t57 = sin(qJ(1));
t58 = cos(qJ(1));
t61 = t58 * pkin(1) + t57 * qJ(2);
t53 = pkin(8) + qJ(4);
t55 = cos(pkin(8));
t60 = t55 * pkin(3) + t63;
t50 = t57 * pkin(1);
t59 = -t58 * qJ(2) + t50;
t41 = g(1) * t57 - g(2) * t58;
t52 = -pkin(7) + t56;
t47 = qJ(5) + t53;
t46 = cos(t53);
t45 = sin(t53);
t44 = cos(t47);
t43 = sin(t47);
t42 = g(1) * t58 + g(2) * t57;
t40 = pkin(4) * t45 + t62;
t1 = [0, 0, 0, 0, 0, 0, -t42, t41, -g(3), -t64, 0, 0, 0, 0, 0, 0, -g(3), t42, -t41, -g(1) * t61 - g(2) * t59 - t64, 0, 0, 0, 0, 0, 0, -g(3) * t55 - t41 * t54, g(3) * t54 - t41 * t55, -t42, -g(1) * (t58 * qJ(3) + t61) - g(2) * (t57 * qJ(3) + t59) - g(3) * t63, 0, 0, 0, 0, 0, 0, -g(3) * t46 - t41 * t45, g(3) * t45 - t41 * t46, -t42, -g(1) * (-t58 * t56 + t57 * t62 + t61) - g(2) * (-t57 * t56 + t50 + (-qJ(2) - t62) * t58) - g(3) * t60, 0, 0, 0, 0, 0, 0, -g(3) * t44 - t41 * t43, g(3) * t43 - t41 * t44, -t42, -g(1) * (t57 * t40 - t58 * t52 + t61) - g(2) * (-t57 * t52 + t50 + (-qJ(2) - t40) * t58) - g(3) * (pkin(4) * t46 + t60);];
U_reg = t1;

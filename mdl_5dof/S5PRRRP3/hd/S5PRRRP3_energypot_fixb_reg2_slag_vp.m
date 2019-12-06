% Calculate inertial parameters regressor of potential energy for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:12
% EndTime: 2019-12-05 16:44:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (110->43), mult. (86->45), div. (0->0), fcn. (73->8), ass. (0->25)
t61 = -pkin(7) - pkin(6);
t58 = pkin(5) + qJ(1);
t65 = g(3) * t58;
t60 = cos(qJ(3));
t44 = t60 * pkin(3) + pkin(2);
t59 = sin(qJ(3));
t64 = t59 * pkin(3) + t58;
t54 = pkin(8) + qJ(2);
t45 = sin(t54);
t46 = cos(t54);
t63 = g(1) * t46 + g(2) * t45;
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t62 = -g(1) * t57 - g(2) * t56;
t55 = qJ(3) + qJ(4);
t53 = -qJ(5) + t61;
t50 = t57 * pkin(1);
t49 = t56 * pkin(1);
t48 = cos(t55);
t47 = sin(t55);
t42 = pkin(4) * t48 + t44;
t41 = g(1) * t45 - g(2) * t46;
t40 = -g(3) * t47 - t63 * t48;
t39 = -g(3) * t48 + t63 * t47;
t1 = [0, 0, 0, 0, 0, 0, t62, g(1) * t56 - g(2) * t57, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t63, t41, -g(3), t62 * pkin(1) - t65, 0, 0, 0, 0, 0, 0, -g(3) * t59 - t63 * t60, -g(3) * t60 + t63 * t59, -t41, -g(1) * (t46 * pkin(2) + t45 * pkin(6) + t50) - g(2) * (t45 * pkin(2) - t46 * pkin(6) + t49) - t65, 0, 0, 0, 0, 0, 0, t40, t39, -t41, -g(1) * (t46 * t44 - t45 * t61 + t50) - g(2) * (t45 * t44 + t46 * t61 + t49) - g(3) * t64, 0, 0, 0, 0, 0, 0, t40, t39, -t41, -g(1) * (t46 * t42 - t45 * t53 + t50) - g(2) * (t45 * t42 + t46 * t53 + t49) - g(3) * (pkin(4) * t47 + t64);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RPRPR10
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:04
% EndTime: 2019-12-31 18:26:04
% DurationCPUTime: 0.07s
% Computational Cost: add. (101->42), mult. (115->51), div. (0->0), fcn. (120->8), ass. (0->28)
t64 = g(3) * pkin(5);
t63 = -pkin(6) + pkin(5);
t62 = g(3) * (-qJ(4) + t63);
t47 = sin(qJ(3));
t48 = sin(qJ(1));
t61 = t48 * t47;
t51 = cos(qJ(1));
t60 = t51 * pkin(1) + t48 * qJ(2);
t59 = qJ(3) + pkin(8);
t43 = t48 * pkin(1);
t58 = -t51 * qJ(2) + t43;
t50 = cos(qJ(3));
t41 = t50 * pkin(3) + pkin(2);
t57 = pkin(3) * t61 + t51 * t41 + t60;
t56 = cos(t59);
t55 = sin(t59);
t29 = -t48 * t55 - t51 * t56;
t30 = -t48 * t56 + t51 * t55;
t54 = g(1) * t30 - g(2) * t29;
t53 = g(1) * t29 + g(2) * t30;
t52 = t48 * t41 + t43 + (-pkin(3) * t47 - qJ(2)) * t51;
t49 = cos(qJ(5));
t46 = sin(qJ(5));
t34 = -g(1) * t51 - g(2) * t48;
t33 = g(1) * t48 - g(2) * t51;
t32 = -t51 * t47 + t48 * t50;
t31 = -t51 * t50 - t61;
t1 = [0, 0, 0, 0, 0, 0, t34, t33, -g(3), -t64, 0, 0, 0, 0, 0, 0, t34, -g(3), -t33, -g(1) * t60 - g(2) * t58 - t64, 0, 0, 0, 0, 0, 0, g(1) * t31 - g(2) * t32, -g(1) * t32 - g(2) * t31, g(3), -g(1) * (t51 * pkin(2) + t60) - g(2) * (t48 * pkin(2) + t58) - g(3) * t63, 0, 0, 0, 0, 0, 0, t53, t54, g(3), -g(1) * t57 - g(2) * t52 - t62, 0, 0, 0, 0, 0, 0, g(3) * t46 + t53 * t49, g(3) * t49 - t53 * t46, -t54, -g(1) * (-t29 * pkin(4) + t30 * pkin(7) + t57) - g(2) * (-t30 * pkin(4) - t29 * pkin(7) + t52) - t62;];
U_reg = t1;

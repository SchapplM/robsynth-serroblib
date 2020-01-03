% Calculate inertial parameters regressor of potential energy for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:23
% EndTime: 2019-12-31 17:32:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->41), mult. (130->47), div. (0->0), fcn. (141->8), ass. (0->25)
t51 = -pkin(5) + qJ(1);
t62 = g(3) * t51;
t61 = cos(qJ(3));
t60 = sin(qJ(3));
t59 = g(3) * qJ(1);
t49 = cos(pkin(7));
t57 = sin(pkin(7));
t58 = t49 * pkin(1) + t57 * qJ(2);
t56 = t49 * pkin(2) + t58;
t55 = t57 * pkin(1) - t49 * qJ(2);
t54 = t57 * pkin(2) + t55;
t31 = -t49 * t61 - t57 * t60;
t32 = t49 * t60 - t57 * t61;
t53 = g(1) * t32 - g(2) * t31;
t52 = g(1) * t31 + g(2) * t32;
t50 = -pkin(6) - qJ(4);
t48 = cos(pkin(8));
t47 = sin(pkin(8));
t46 = pkin(8) + qJ(5);
t40 = cos(t46);
t39 = sin(t46);
t38 = t48 * pkin(4) + pkin(3);
t34 = -g(1) * t49 - g(2) * t57;
t33 = g(1) * t57 - g(2) * t49;
t1 = [0, 0, 0, 0, 0, 0, t34, t33, -g(3), -t59, 0, 0, 0, 0, 0, 0, t34, -g(3), -t33, -g(1) * t58 - g(2) * t55 - t59, 0, 0, 0, 0, 0, 0, t52, t53, g(3), -g(1) * t56 - g(2) * t54 - t62, 0, 0, 0, 0, 0, 0, g(3) * t47 + t52 * t48, g(3) * t48 - t52 * t47, -t53, -g(1) * (-t31 * pkin(3) + t32 * qJ(4) + t56) - g(2) * (-t32 * pkin(3) - t31 * qJ(4) + t54) - t62, 0, 0, 0, 0, 0, 0, g(3) * t39 + t52 * t40, g(3) * t40 - t52 * t39, -t53, -g(1) * (-t31 * t38 - t32 * t50 + t56) - g(2) * (t31 * t50 - t32 * t38 + t54) - g(3) * (-t47 * pkin(4) + t51);];
U_reg = t1;

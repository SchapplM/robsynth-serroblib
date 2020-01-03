% Calculate inertial parameters regressor of potential energy for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (99->40), mult. (78->41), div. (0->0), fcn. (65->8), ass. (0->22)
t47 = sin(pkin(8));
t59 = pkin(4) * t47;
t49 = qJ(2) + pkin(5);
t58 = g(3) * t49;
t46 = qJ(1) + pkin(7);
t40 = sin(t46);
t51 = sin(qJ(1));
t57 = t51 * pkin(1) + t40 * pkin(2);
t56 = pkin(3) + t49;
t42 = cos(t46);
t52 = cos(qJ(1));
t55 = t52 * pkin(1) + t42 * pkin(2) + t40 * qJ(3);
t54 = -t42 * qJ(3) + t57;
t34 = g(1) * t40 - g(2) * t42;
t53 = -g(1) * t52 - g(2) * t51;
t50 = -pkin(6) - qJ(4);
t48 = cos(pkin(8));
t45 = pkin(8) + qJ(5);
t41 = cos(t45);
t39 = sin(t45);
t35 = g(1) * t42 + g(2) * t40;
t1 = [0, 0, 0, 0, 0, 0, t53, g(1) * t51 - g(2) * t52, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t35, t34, -g(3), t53 * pkin(1) - t58, 0, 0, 0, 0, 0, 0, -g(3), t35, -t34, -g(1) * t55 - g(2) * t54 - t58, 0, 0, 0, 0, 0, 0, -g(3) * t48 - t34 * t47, g(3) * t47 - t34 * t48, -t35, -g(1) * (t42 * qJ(4) + t55) - g(2) * (t40 * qJ(4) + t54) - g(3) * t56, 0, 0, 0, 0, 0, 0, -g(3) * t41 - t34 * t39, g(3) * t39 - t34 * t41, -t35, -g(1) * (t40 * t59 - t42 * t50 + t55) - g(2) * (-t40 * t50 + (-qJ(3) - t59) * t42 + t57) - g(3) * (t48 * pkin(4) + t56);];
U_reg = t1;

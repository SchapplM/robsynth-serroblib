% Calculate inertial parameters regressor of potential energy for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:40
% EndTime: 2019-12-31 16:52:40
% DurationCPUTime: 0.05s
% Computational Cost: add. (67->34), mult. (71->41), div. (0->0), fcn. (61->8), ass. (0->19)
t58 = g(3) * pkin(4);
t51 = sin(pkin(7));
t57 = t51 * pkin(2) + pkin(4);
t52 = cos(pkin(7));
t42 = t52 * pkin(2) + pkin(1);
t53 = -pkin(5) - qJ(2);
t50 = pkin(7) + qJ(3);
t54 = sin(qJ(1));
t55 = cos(qJ(1));
t56 = g(1) * t55 + g(2) * t54;
t49 = -pkin(6) + t53;
t45 = qJ(4) + t50;
t44 = cos(t50);
t43 = sin(t50);
t41 = cos(t45);
t40 = sin(t45);
t39 = g(1) * t54 - g(2) * t55;
t38 = pkin(3) * t44 + t42;
t1 = [0, 0, 0, 0, 0, 0, -t56, t39, -g(3), -t58, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t56 * t52, -g(3) * t52 + t56 * t51, -t39, -g(1) * (t55 * pkin(1) + t54 * qJ(2)) - g(2) * (t54 * pkin(1) - t55 * qJ(2)) - t58, 0, 0, 0, 0, 0, 0, -g(3) * t43 - t56 * t44, -g(3) * t44 + t56 * t43, -t39, -g(1) * (t55 * t42 - t54 * t53) - g(2) * (t54 * t42 + t55 * t53) - g(3) * t57, 0, 0, 0, 0, 0, 0, -g(3) * t40 - t56 * t41, -g(3) * t41 + t56 * t40, -t39, -g(1) * (t55 * t38 - t54 * t49) - g(2) * (t54 * t38 + t55 * t49) - g(3) * (pkin(3) * t43 + t57);];
U_reg = t1;

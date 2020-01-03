% Calculate inertial parameters regressor of potential energy for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:34
% EndTime: 2019-12-31 16:50:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (77->36), mult. (79->48), div. (0->0), fcn. (73->8), ass. (0->22)
t42 = qJ(2) + pkin(4);
t57 = g(3) * t42;
t44 = sin(qJ(3));
t56 = g(3) * t44;
t43 = sin(qJ(4));
t47 = cos(qJ(3));
t55 = t43 * t47;
t46 = cos(qJ(4));
t54 = t46 * t47;
t41 = qJ(1) + pkin(7);
t37 = sin(t41);
t38 = cos(t41);
t48 = cos(qJ(1));
t53 = t48 * pkin(1) + t38 * pkin(2) + t37 * pkin(5);
t45 = sin(qJ(1));
t52 = t45 * pkin(1) + t37 * pkin(2) - t38 * pkin(5);
t51 = pkin(3) * t47 + pkin(6) * t44;
t50 = g(1) * t38 + g(2) * t37;
t49 = -g(1) * t48 - g(2) * t45;
t33 = g(1) * t37 - g(2) * t38;
t32 = -g(3) * t47 + t50 * t44;
t1 = [0, 0, 0, 0, 0, 0, t49, g(1) * t45 - g(2) * t48, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t50, t33, -g(3), t49 * pkin(1) - t57, 0, 0, 0, 0, 0, 0, -t50 * t47 - t56, t32, -t33, -g(1) * t53 - g(2) * t52 - t57, 0, 0, 0, 0, 0, 0, -g(1) * (t37 * t43 + t38 * t54) - g(2) * (t37 * t54 - t38 * t43) - t46 * t56, -g(1) * (t37 * t46 - t38 * t55) - g(2) * (-t37 * t55 - t38 * t46) + t43 * t56, -t32, -g(1) * (t51 * t38 + t53) - g(2) * (t51 * t37 + t52) - g(3) * (t44 * pkin(3) - t47 * pkin(6) + t42);];
U_reg = t1;

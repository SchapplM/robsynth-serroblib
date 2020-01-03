% Calculate inertial parameters regressor of potential energy for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (65->32), mult. (78->37), div. (0->0), fcn. (68->6), ass. (0->19)
t56 = g(3) * pkin(4);
t46 = sin(pkin(6));
t55 = t46 * pkin(2) + pkin(4);
t47 = cos(pkin(6));
t39 = t47 * pkin(2) + pkin(1);
t48 = -pkin(5) - qJ(2);
t49 = sin(qJ(1));
t50 = cos(qJ(1));
t54 = t49 * t39 + t50 * t48;
t53 = t50 * t39 - t49 * t48;
t52 = g(1) * t50 + g(2) * t49;
t45 = pkin(6) + qJ(3);
t41 = sin(t45);
t42 = cos(t45);
t51 = pkin(3) * t42 + qJ(4) * t41;
t38 = g(1) * t49 - g(2) * t50;
t35 = -g(3) * t41 - t52 * t42;
t34 = -g(3) * t42 + t52 * t41;
t1 = [0, 0, 0, 0, 0, 0, -t52, t38, -g(3), -t56, 0, 0, 0, 0, 0, 0, -g(3) * t46 - t52 * t47, -g(3) * t47 + t52 * t46, -t38, -g(1) * (t50 * pkin(1) + t49 * qJ(2)) - g(2) * (t49 * pkin(1) - t50 * qJ(2)) - t56, 0, 0, 0, 0, 0, 0, t35, t34, -t38, -g(1) * t53 - g(2) * t54 - g(3) * t55, 0, 0, 0, 0, 0, 0, t35, -t38, -t34, -g(1) * (t51 * t50 + t53) - g(2) * (t51 * t49 + t54) - g(3) * (t41 * pkin(3) - t42 * qJ(4) + t55);];
U_reg = t1;

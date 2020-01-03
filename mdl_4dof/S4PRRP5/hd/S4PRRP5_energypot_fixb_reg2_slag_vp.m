% Calculate inertial parameters regressor of potential energy for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:21
% EndTime: 2019-12-31 16:29:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (61->41), mult. (115->54), div. (0->0), fcn. (113->6), ass. (0->24)
t43 = cos(qJ(3));
t33 = t43 * pkin(3) + pkin(2);
t40 = -qJ(4) - pkin(5);
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t56 = t33 * t44 - t40 * t42;
t55 = g(3) * t42;
t54 = g(3) * qJ(1);
t38 = sin(pkin(6));
t41 = sin(qJ(3));
t52 = t38 * t41;
t50 = t41 * t44;
t49 = t43 * t44;
t39 = cos(pkin(6));
t48 = t39 * pkin(1) + t38 * pkin(4);
t35 = t38 * pkin(1);
t47 = -t39 * pkin(4) + t35;
t46 = pkin(2) * t44 + pkin(5) * t42;
t45 = g(1) * t39 + g(2) * t38;
t32 = g(1) * t38 - g(2) * t39;
t31 = -g(3) * t44 + t45 * t42;
t30 = -g(1) * (t39 * t49 + t52) - g(2) * (t38 * t49 - t39 * t41) - t43 * t55;
t29 = -g(1) * (t38 * t43 - t39 * t50) - g(2) * (-t38 * t50 - t39 * t43) + t41 * t55;
t1 = [0, 0, 0, 0, 0, 0, -t45, t32, -g(3), -t54, 0, 0, 0, 0, 0, 0, -t45 * t44 - t55, t31, -t32, -g(1) * t48 - g(2) * t47 - t54, 0, 0, 0, 0, 0, 0, t30, t29, -t31, -g(1) * (t46 * t39 + t48) - g(2) * (t46 * t38 + t47) - g(3) * (t42 * pkin(2) - t44 * pkin(5) + qJ(1)), 0, 0, 0, 0, 0, 0, t30, t29, -t31, -g(1) * (pkin(3) * t52 + t48) - g(2) * (t56 * t38 + t35) - g(3) * (t42 * t33 + t44 * t40 + qJ(1)) + (-g(1) * t56 - g(2) * (-pkin(3) * t41 - pkin(4))) * t39;];
U_reg = t1;

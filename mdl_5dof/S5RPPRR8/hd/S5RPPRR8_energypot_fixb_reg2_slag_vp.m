% Calculate inertial parameters regressor of potential energy for
% S5RPPRR8
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:13
% EndTime: 2019-12-31 18:01:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (101->42), mult. (115->51), div. (0->0), fcn. (120->8), ass. (0->28)
t65 = g(3) * pkin(5);
t62 = -qJ(3) + pkin(5);
t64 = g(3) * (-pkin(6) + t62);
t47 = sin(pkin(8));
t50 = sin(qJ(1));
t63 = t50 * t47;
t52 = cos(qJ(1));
t61 = t52 * pkin(1) + t50 * qJ(2);
t60 = pkin(8) + qJ(4);
t44 = t50 * pkin(1);
t59 = -t52 * qJ(2) + t44;
t48 = cos(pkin(8));
t42 = t48 * pkin(3) + pkin(2);
t58 = pkin(3) * t63 + t52 * t42 + t61;
t57 = cos(t60);
t56 = sin(t60);
t30 = -t50 * t56 - t52 * t57;
t31 = -t50 * t57 + t52 * t56;
t55 = g(1) * t31 - g(2) * t30;
t54 = g(1) * t30 + g(2) * t31;
t53 = t50 * t42 + t44 + (-pkin(3) * t47 - qJ(2)) * t52;
t51 = cos(qJ(5));
t49 = sin(qJ(5));
t37 = -g(1) * t52 - g(2) * t50;
t36 = g(1) * t50 - g(2) * t52;
t33 = -t52 * t47 + t50 * t48;
t32 = -t52 * t48 - t63;
t1 = [0, 0, 0, 0, 0, 0, t37, t36, -g(3), -t65, 0, 0, 0, 0, 0, 0, t37, -g(3), -t36, -g(1) * t61 - g(2) * t59 - t65, 0, 0, 0, 0, 0, 0, g(1) * t32 - g(2) * t33, -g(1) * t33 - g(2) * t32, g(3), -g(1) * (t52 * pkin(2) + t61) - g(2) * (t50 * pkin(2) + t59) - g(3) * t62, 0, 0, 0, 0, 0, 0, t54, t55, g(3), -g(1) * t58 - g(2) * t53 - t64, 0, 0, 0, 0, 0, 0, g(3) * t49 + t54 * t51, g(3) * t51 - t54 * t49, -t55, -g(1) * (-t30 * pkin(4) + t31 * pkin(7) + t58) - g(2) * (-t31 * pkin(4) - t30 * pkin(7) + t53) - t64;];
U_reg = t1;

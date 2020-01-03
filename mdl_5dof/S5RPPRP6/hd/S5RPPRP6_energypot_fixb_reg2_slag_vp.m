% Calculate inertial parameters regressor of potential energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:16
% DurationCPUTime: 0.07s
% Computational Cost: add. (81->43), mult. (99->44), div. (0->0), fcn. (86->6), ass. (0->24)
t49 = pkin(7) + qJ(4);
t43 = sin(t49);
t44 = cos(t49);
t66 = pkin(4) * t43 - qJ(5) * t44;
t65 = g(3) * pkin(5);
t64 = pkin(2) + pkin(5);
t50 = sin(pkin(7));
t63 = pkin(3) * t50;
t53 = sin(qJ(1));
t54 = cos(qJ(1));
t61 = t54 * pkin(1) + t53 * qJ(2);
t51 = cos(pkin(7));
t59 = t51 * pkin(3) + t64;
t58 = t53 * t63 + t61;
t57 = -qJ(2) - t63;
t47 = t53 * pkin(1);
t52 = -pkin(6) - qJ(3);
t56 = -t53 * t52 + t47;
t55 = -t54 * qJ(2) + t47;
t40 = g(1) * t53 - g(2) * t54;
t41 = g(1) * t54 + g(2) * t53;
t39 = -g(3) * t43 + t40 * t44;
t38 = -g(3) * t44 - t40 * t43;
t1 = [0, 0, 0, 0, 0, 0, -t41, t40, -g(3), -t65, 0, 0, 0, 0, 0, 0, -g(3), t41, -t40, -g(1) * t61 - g(2) * t55 - t65, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t40 * t50, g(3) * t50 - t40 * t51, -t41, -g(1) * (t54 * qJ(3) + t61) - g(2) * (t53 * qJ(3) + t55) - g(3) * t64, 0, 0, 0, 0, 0, 0, t38, -t39, -t41, -g(1) * (-t54 * t52 + t58) - g(2) * (t57 * t54 + t56) - g(3) * t59, 0, 0, 0, 0, 0, 0, t38, -t41, t39, -g(1) * (t66 * t53 + t58) - g(2) * t56 - g(3) * (t44 * pkin(4) + t43 * qJ(5) + t59) + (g(1) * t52 - g(2) * (t57 - t66)) * t54;];
U_reg = t1;

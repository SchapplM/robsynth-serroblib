% Calculate inertial parameters regressor of potential energy for
% S5RPPRR12
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:20
% EndTime: 2019-12-31 18:07:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (89->50), mult. (112->58), div. (0->0), fcn. (103->8), ass. (0->30)
t50 = pkin(8) + qJ(4);
t44 = sin(t50);
t45 = cos(t50);
t74 = pkin(4) * t44 - pkin(7) * t45;
t73 = g(3) * pkin(5);
t72 = pkin(2) + pkin(5);
t51 = sin(pkin(8));
t71 = pkin(3) * t51;
t68 = g(3) * t45;
t54 = sin(qJ(5));
t55 = sin(qJ(1));
t67 = t55 * t54;
t56 = cos(qJ(5));
t66 = t55 * t56;
t57 = cos(qJ(1));
t65 = t57 * t54;
t64 = t57 * t56;
t63 = t57 * pkin(1) + t55 * qJ(2);
t52 = cos(pkin(8));
t62 = t52 * pkin(3) + t72;
t61 = t55 * t71 + t63;
t60 = -qJ(2) - t71;
t48 = t55 * pkin(1);
t53 = -pkin(6) - qJ(3);
t59 = -t55 * t53 + t48;
t58 = -t57 * qJ(2) + t48;
t41 = g(1) * t55 - g(2) * t57;
t42 = g(1) * t57 + g(2) * t55;
t40 = -g(3) * t44 + t41 * t45;
t1 = [0, 0, 0, 0, 0, 0, -t42, t41, -g(3), -t73, 0, 0, 0, 0, 0, 0, -g(3), t42, -t41, -g(1) * t63 - g(2) * t58 - t73, 0, 0, 0, 0, 0, 0, -g(3) * t52 - t41 * t51, g(3) * t51 - t41 * t52, -t42, -g(1) * (t57 * qJ(3) + t63) - g(2) * (t55 * qJ(3) + t58) - g(3) * t72, 0, 0, 0, 0, 0, 0, -t41 * t44 - t68, -t40, -t42, -g(1) * (-t57 * t53 + t61) - g(2) * (t60 * t57 + t59) - g(3) * t62, 0, 0, 0, 0, 0, 0, -g(1) * (t44 * t66 + t65) - g(2) * (-t44 * t64 + t67) - t56 * t68, -g(1) * (-t44 * t67 + t64) - g(2) * (t44 * t65 + t66) + t54 * t68, t40, -g(1) * (t55 * t74 + t61) - g(2) * t59 - g(3) * (t45 * pkin(4) + t44 * pkin(7) + t62) + (g(1) * t53 - g(2) * (t60 - t74)) * t57;];
U_reg = t1;

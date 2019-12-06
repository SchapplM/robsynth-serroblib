% Calculate inertial parameters regressor of potential energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:18
% EndTime: 2019-12-05 15:41:18
% DurationCPUTime: 0.09s
% Computational Cost: add. (88->53), mult. (171->61), div. (0->0), fcn. (170->6), ass. (0->34)
t54 = sin(pkin(7));
t57 = sin(qJ(2));
t68 = qJ(3) * t57;
t59 = cos(qJ(2));
t73 = t54 * t59;
t76 = pkin(2) * t73 + t54 * t68;
t75 = g(3) * t59;
t74 = g(3) * qJ(1);
t55 = cos(pkin(7));
t72 = t55 * t59;
t56 = sin(qJ(4));
t71 = t56 * t57;
t58 = cos(qJ(4));
t70 = t57 * t58;
t69 = t55 * pkin(1) + t54 * pkin(5);
t67 = t57 * pkin(2) + qJ(1);
t49 = t54 * pkin(1);
t66 = -t55 * pkin(5) + t49;
t65 = pkin(2) * t72 + t55 * t68 + t69;
t64 = g(1) * t55 + g(2) * t54;
t63 = -t59 * qJ(3) + t67;
t62 = t54 * pkin(3) + pkin(6) * t72 + t65;
t61 = pkin(6) * t73 + t49 + (-pkin(3) - pkin(5)) * t55 + t76;
t35 = t54 * t56 - t55 * t70;
t37 = t54 * t70 + t55 * t56;
t60 = g(1) * t35 - g(2) * t37 + t58 * t75;
t52 = t57 * pkin(6);
t39 = g(1) * t54 - g(2) * t55;
t38 = t54 * t71 - t55 * t58;
t36 = t54 * t58 + t55 * t71;
t34 = g(3) * t57 + t64 * t59;
t33 = t64 * t57 - t75;
t32 = -g(1) * t36 - g(2) * t38 + t56 * t75;
t1 = [0, 0, 0, 0, 0, 0, -t64, t39, -g(3), -t74, 0, 0, 0, 0, 0, 0, -t34, t33, -t39, -g(1) * t69 - g(2) * t66 - t74, 0, 0, 0, 0, 0, 0, -t39, t34, -t33, -g(1) * t65 - g(2) * (t66 + t76) - g(3) * t63, 0, 0, 0, 0, 0, 0, t32, t60, -t34, -g(1) * t62 - g(2) * t61 - g(3) * (t52 + t63), 0, 0, 0, 0, 0, 0, t32, -t34, -t60, -g(1) * (t36 * pkin(4) + t35 * qJ(5) + t62) - g(2) * (t38 * pkin(4) - t37 * qJ(5) + t61) - g(3) * (t52 + t67) - (-pkin(4) * t56 + qJ(5) * t58 - qJ(3)) * t75;];
U_reg = t1;

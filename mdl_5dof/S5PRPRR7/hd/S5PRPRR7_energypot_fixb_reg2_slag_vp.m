% Calculate inertial parameters regressor of potential energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:43
% EndTime: 2019-12-05 16:00:44
% DurationCPUTime: 0.16s
% Computational Cost: add. (96->65), mult. (161->82), div. (0->0), fcn. (156->8), ass. (0->35)
t50 = sin(pkin(8));
t53 = sin(qJ(2));
t64 = qJ(3) * t53;
t55 = cos(qJ(2));
t71 = t50 * t55;
t75 = pkin(2) * t71 + t50 * t64;
t74 = g(3) * t55;
t73 = g(3) * qJ(1);
t72 = t50 * t53;
t51 = cos(pkin(8));
t70 = t51 * t53;
t69 = t51 * t55;
t52 = sin(qJ(4));
t68 = t52 * t53;
t54 = cos(qJ(4));
t67 = t53 * t54;
t56 = -pkin(7) - pkin(6);
t66 = t55 * t56;
t65 = t51 * pkin(1) + t50 * pkin(5);
t63 = t53 * pkin(2) + qJ(1);
t62 = t50 * t68;
t45 = t50 * pkin(1);
t61 = t45 + t75;
t60 = -t51 * pkin(5) + t45;
t59 = pkin(2) * t69 + t51 * t64 + t65;
t58 = g(1) * t51 + g(2) * t50;
t57 = -t55 * qJ(3) + t63;
t49 = qJ(4) + qJ(5);
t43 = cos(t49);
t42 = sin(t49);
t41 = t54 * pkin(4) + pkin(3);
t36 = g(1) * t50 - g(2) * t51;
t35 = g(3) * t53 + t58 * t55;
t34 = t58 * t53 - t74;
t1 = [0, 0, 0, 0, 0, 0, -t58, t36, -g(3), -t73, 0, 0, 0, 0, 0, 0, -t35, t34, -t36, -g(1) * t65 - g(2) * t60 - t73, 0, 0, 0, 0, 0, 0, -t36, t35, -t34, -g(1) * t59 - g(2) * (t60 + t75) - g(3) * t57, 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t54 + t51 * t68) - g(2) * (-t51 * t54 + t62) + t52 * t74, -g(1) * (-t50 * t52 + t51 * t67) - g(2) * (t50 * t67 + t51 * t52) + t54 * t74, -t35, -g(1) * (t50 * pkin(3) + pkin(6) * t69 + t59) - g(2) * (pkin(6) * t71 + (-pkin(3) - pkin(5)) * t51 + t61) - g(3) * (t53 * pkin(6) + t57), 0, 0, 0, 0, 0, 0, -g(1) * (t42 * t70 + t50 * t43) - g(2) * (t42 * t72 - t51 * t43) + t42 * t74, -g(1) * (-t50 * t42 + t43 * t70) - g(2) * (t51 * t42 + t43 * t72) + t43 * t74, -t35, -g(1) * (t50 * t41 + t59) - g(2) * (pkin(4) * t62 - t50 * t66 + t61) - g(3) * (-t53 * t56 + (-pkin(4) * t52 - qJ(3)) * t55 + t63) + (-g(1) * (pkin(4) * t68 - t66) - g(2) * (-pkin(5) - t41)) * t51;];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:06
% EndTime: 2019-12-31 17:58:07
% DurationCPUTime: 0.08s
% Computational Cost: add. (130->49), mult. (106->59), div. (0->0), fcn. (97->10), ass. (0->32)
t58 = pkin(9) + qJ(4);
t51 = sin(t58);
t79 = g(3) * t51;
t62 = qJ(2) + pkin(5);
t78 = g(3) * t62;
t59 = qJ(1) + pkin(8);
t52 = sin(t59);
t64 = sin(qJ(5));
t77 = t52 * t64;
t66 = cos(qJ(5));
t76 = t52 * t66;
t54 = cos(t59);
t75 = t54 * t64;
t74 = t54 * t66;
t61 = cos(pkin(9));
t50 = t61 * pkin(3) + pkin(2);
t65 = sin(qJ(1));
t56 = t65 * pkin(1);
t63 = -pkin(6) - qJ(3);
t73 = t52 * t50 + t54 * t63 + t56;
t60 = sin(pkin(9));
t72 = t60 * pkin(3) + t62;
t67 = cos(qJ(1));
t57 = t67 * pkin(1);
t71 = t54 * t50 - t52 * t63 + t57;
t53 = cos(t58);
t70 = pkin(4) * t53 + pkin(7) * t51;
t69 = g(1) * t54 + g(2) * t52;
t68 = -g(1) * t67 - g(2) * t65;
t45 = g(1) * t52 - g(2) * t54;
t44 = -g(3) * t53 + t69 * t51;
t1 = [0, 0, 0, 0, 0, 0, t68, g(1) * t65 - g(2) * t67, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t69, t45, -g(3), t68 * pkin(1) - t78, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t69 * t61, -g(3) * t61 + t69 * t60, -t45, -g(1) * (t54 * pkin(2) + t52 * qJ(3) + t57) - g(2) * (t52 * pkin(2) - t54 * qJ(3) + t56) - t78, 0, 0, 0, 0, 0, 0, -t69 * t53 - t79, t44, -t45, -g(1) * t71 - g(2) * t73 - g(3) * t72, 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t74 + t77) - g(2) * (t53 * t76 - t75) - t66 * t79, -g(1) * (-t53 * t75 + t76) - g(2) * (-t53 * t77 - t74) + t64 * t79, -t44, -g(1) * (t70 * t54 + t71) - g(2) * (t70 * t52 + t73) - g(3) * (t51 * pkin(4) - t53 * pkin(7) + t72);];
U_reg = t1;

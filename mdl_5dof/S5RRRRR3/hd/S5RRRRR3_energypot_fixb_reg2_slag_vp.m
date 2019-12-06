% Calculate inertial parameters regressor of potential energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:28
% EndTime: 2019-12-05 18:56:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->49), mult. (131->66), div. (0->0), fcn. (128->10), ass. (0->35)
t59 = qJ(2) + qJ(3);
t53 = sin(t59);
t64 = cos(qJ(2));
t86 = pkin(1) * t64 + pkin(5) * t53;
t85 = g(3) * pkin(4);
t82 = g(3) * t53;
t61 = sin(qJ(2));
t81 = t61 * pkin(1) + pkin(4);
t55 = cos(t59);
t62 = sin(qJ(1));
t80 = t55 * t62;
t65 = cos(qJ(1));
t79 = t55 * t65;
t58 = qJ(4) + qJ(5);
t52 = sin(t58);
t78 = t62 * t52;
t54 = cos(t58);
t77 = t62 * t54;
t60 = sin(qJ(4));
t76 = t62 * t60;
t63 = cos(qJ(4));
t75 = t62 * t63;
t74 = t65 * t52;
t73 = t65 * t54;
t72 = t65 * t60;
t71 = t65 * t63;
t70 = t86 * t62;
t69 = t86 * t65;
t68 = -t55 * pkin(5) + t81;
t67 = g(1) * t65 + g(2) * t62;
t66 = t67 * t64;
t51 = t63 * pkin(3) + pkin(2);
t45 = g(1) * t62 - g(2) * t65;
t44 = -g(3) * t55 + t67 * t53;
t1 = [0, 0, 0, 0, 0, 0, -t67, t45, -g(3), -t85, 0, 0, 0, 0, 0, 0, -g(3) * t61 - t66, -g(3) * t64 + t67 * t61, -t45, -t85, 0, 0, 0, 0, 0, 0, -t67 * t55 - t82, t44, -t45, -pkin(1) * t66 - g(3) * t81, 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t71 + t76) - g(2) * (t55 * t75 - t72) - t63 * t82, -g(1) * (-t55 * t72 + t75) - g(2) * (-t55 * t76 - t71) + t60 * t82, -t44, -g(1) * (pkin(2) * t79 + t69) - g(2) * (pkin(2) * t80 + t70) - g(3) * (t53 * pkin(2) + t68), 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t73 + t78) - g(2) * (t55 * t77 - t74) - t54 * t82, -g(1) * (-t55 * t74 + t77) - g(2) * (-t55 * t78 - t73) + t52 * t82, -t44, -g(1) * (pkin(3) * t76 + t51 * t79 + t69) - g(2) * (-pkin(3) * t72 + t51 * t80 + t70) - g(3) * (t53 * t51 + t68);];
U_reg = t1;

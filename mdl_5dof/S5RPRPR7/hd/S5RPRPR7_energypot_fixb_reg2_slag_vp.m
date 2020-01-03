% Calculate inertial parameters regressor of potential energy for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:50
% EndTime: 2019-12-31 18:19:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (130->49), mult. (106->59), div. (0->0), fcn. (97->10), ass. (0->32)
t61 = qJ(3) + pkin(9);
t54 = sin(t61);
t82 = g(3) * t54;
t64 = qJ(2) + pkin(5);
t81 = g(3) * t64;
t62 = qJ(1) + pkin(8);
t55 = sin(t62);
t65 = sin(qJ(5));
t80 = t55 * t65;
t68 = cos(qJ(5));
t79 = t55 * t68;
t57 = cos(t62);
t78 = t57 * t65;
t77 = t57 * t68;
t69 = cos(qJ(3));
t53 = pkin(3) * t69 + pkin(2);
t67 = sin(qJ(1));
t59 = t67 * pkin(1);
t63 = -qJ(4) - pkin(6);
t76 = t55 * t53 + t57 * t63 + t59;
t66 = sin(qJ(3));
t75 = t66 * pkin(3) + t64;
t70 = cos(qJ(1));
t60 = t70 * pkin(1);
t74 = t57 * t53 - t55 * t63 + t60;
t56 = cos(t61);
t73 = pkin(4) * t56 + pkin(7) * t54;
t72 = g(1) * t57 + g(2) * t55;
t71 = -g(1) * t70 - g(2) * t67;
t48 = g(1) * t55 - g(2) * t57;
t47 = -g(3) * t56 + t54 * t72;
t1 = [0, 0, 0, 0, 0, 0, t71, g(1) * t67 - g(2) * t70, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t72, t48, -g(3), pkin(1) * t71 - t81, 0, 0, 0, 0, 0, 0, -g(3) * t66 - t69 * t72, -g(3) * t69 + t66 * t72, -t48, -g(1) * (pkin(2) * t57 + pkin(6) * t55 + t60) - g(2) * (pkin(2) * t55 - pkin(6) * t57 + t59) - t81, 0, 0, 0, 0, 0, 0, -t56 * t72 - t82, t47, -t48, -g(1) * t74 - g(2) * t76 - g(3) * t75, 0, 0, 0, 0, 0, 0, -g(1) * (t56 * t77 + t80) - g(2) * (t56 * t79 - t78) - t68 * t82, -g(1) * (-t56 * t78 + t79) - g(2) * (-t56 * t80 - t77) + t65 * t82, -t47, -g(1) * (t57 * t73 + t74) - g(2) * (t55 * t73 + t76) - g(3) * (pkin(4) * t54 - pkin(7) * t56 + t75);];
U_reg = t1;

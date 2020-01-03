% Calculate inertial parameters regressor of potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:31
% EndTime: 2020-01-03 11:20:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (138->60), mult. (130->82), div. (0->0), fcn. (125->10), ass. (0->32)
t57 = sin(pkin(8));
t77 = g(1) * t57;
t60 = qJ(2) + pkin(5);
t76 = g(1) * t60;
t63 = cos(qJ(1));
t75 = t63 * pkin(1);
t55 = qJ(1) + pkin(7);
t50 = sin(t55);
t59 = cos(pkin(8));
t74 = t50 * t59;
t52 = cos(t55);
t73 = t52 * t59;
t56 = sin(pkin(9));
t72 = t56 * t59;
t58 = cos(pkin(9));
t71 = t58 * t59;
t62 = sin(qJ(1));
t70 = t62 * pkin(1) + t50 * pkin(2);
t69 = qJ(4) * t57;
t68 = pkin(4) * t56 + qJ(3);
t67 = -g(2) * t50 + g(3) * t52;
t66 = -g(2) * t62 + g(3) * t63;
t65 = -t50 * qJ(3) - t75;
t48 = t58 * pkin(4) + pkin(3);
t61 = -pkin(6) - qJ(4);
t64 = t48 * t59 - t57 * t61;
t54 = pkin(9) + qJ(5);
t51 = cos(t54);
t49 = sin(t54);
t46 = g(2) * t52 + g(3) * t50;
t45 = g(1) * t59 + t67 * t57;
t1 = [0, 0, 0, 0, 0, 0, t66, -g(2) * t63 - g(3) * t62, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t67, -t46, -g(1), t66 * pkin(1) - t76, 0, 0, 0, 0, 0, 0, t67 * t59 - t77, -t45, t46, -t76 - g(2) * (-t52 * qJ(3) + t70) - g(3) * (-t52 * pkin(2) + t65), 0, 0, 0, 0, 0, 0, -t58 * t77 - g(2) * (t50 * t71 - t52 * t56) - g(3) * (-t50 * t56 - t52 * t71), t56 * t77 - g(2) * (-t50 * t72 - t52 * t58) - g(3) * (-t50 * t58 + t52 * t72), t45, -g(1) * (t57 * pkin(3) - t59 * qJ(4) + t60) - g(2) * (pkin(3) * t74 + t50 * t69 + t70) - g(3) * t65 + (g(2) * qJ(3) - g(3) * (-pkin(3) * t59 - pkin(2) - t69)) * t52, 0, 0, 0, 0, 0, 0, -t51 * t77 - g(2) * (-t52 * t49 + t51 * t74) - g(3) * (-t50 * t49 - t51 * t73), t49 * t77 - g(2) * (-t49 * t74 - t52 * t51) - g(3) * (t49 * t73 - t50 * t51), t45, -g(1) * (t57 * t48 + t59 * t61 + t60) - g(2) * t70 + g(3) * t75 + (-g(2) * t64 + g(3) * t68) * t50 + (g(2) * t68 - g(3) * (-pkin(2) - t64)) * t52;];
U_reg = t1;

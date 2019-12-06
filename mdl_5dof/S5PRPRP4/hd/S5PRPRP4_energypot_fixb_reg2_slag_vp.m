% Calculate inertial parameters regressor of potential energy for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:05
% EndTime: 2019-12-05 15:36:05
% DurationCPUTime: 0.11s
% Computational Cost: add. (129->49), mult. (156->61), div. (0->0), fcn. (155->8), ass. (0->34)
t58 = qJ(2) + pkin(8);
t54 = sin(t58);
t55 = cos(t58);
t82 = pkin(3) * t55 + pkin(6) * t54;
t79 = g(3) * t54;
t78 = g(3) * qJ(1);
t59 = sin(pkin(7));
t62 = sin(qJ(4));
t77 = t59 * t62;
t64 = cos(qJ(4));
t76 = t59 * t64;
t60 = cos(pkin(7));
t75 = t60 * t62;
t74 = t60 * t64;
t65 = cos(qJ(2));
t53 = t65 * pkin(2) + pkin(1);
t61 = -qJ(3) - pkin(5);
t73 = t59 * t53 + t60 * t61;
t63 = sin(qJ(2));
t72 = t63 * pkin(2) + qJ(1);
t71 = t60 * t53 - t59 * t61;
t70 = t82 * t59 + t73;
t69 = g(1) * t60 + g(2) * t59;
t68 = t54 * pkin(3) - t55 * pkin(6) + t72;
t67 = t82 * t60 + t71;
t38 = t55 * t77 + t74;
t40 = t55 * t75 - t76;
t66 = g(1) * t40 + g(2) * t38 + t62 * t79;
t42 = g(1) * t59 - g(2) * t60;
t41 = t55 * t74 + t77;
t39 = t55 * t76 - t75;
t37 = -g(3) * t55 + t69 * t54;
t36 = -g(1) * t41 - g(2) * t39 - t64 * t79;
t1 = [0, 0, 0, 0, 0, 0, -t69, t42, -g(3), -t78, 0, 0, 0, 0, 0, 0, -g(3) * t63 - t69 * t65, -g(3) * t65 + t69 * t63, -t42, -g(1) * (t60 * pkin(1) + t59 * pkin(5)) - g(2) * (t59 * pkin(1) - t60 * pkin(5)) - t78, 0, 0, 0, 0, 0, 0, -t69 * t55 - t79, t37, -t42, -g(1) * t71 - g(2) * t73 - g(3) * t72, 0, 0, 0, 0, 0, 0, t36, t66, -t37, -g(1) * t67 - g(2) * t70 - g(3) * t68, 0, 0, 0, 0, 0, 0, t36, -t37, -t66, -g(1) * (t41 * pkin(4) + t40 * qJ(5) + t67) - g(2) * (t39 * pkin(4) + t38 * qJ(5) + t70) - g(3) * ((pkin(4) * t64 + qJ(5) * t62) * t54 + t68);];
U_reg = t1;

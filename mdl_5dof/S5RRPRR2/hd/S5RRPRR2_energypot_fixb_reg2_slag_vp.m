% Calculate inertial parameters regressor of potential energy for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:26
% EndTime: 2019-12-05 18:28:26
% DurationCPUTime: 0.09s
% Computational Cost: add. (118->45), mult. (99->54), div. (0->0), fcn. (86->10), ass. (0->25)
t85 = g(3) * pkin(5);
t78 = sin(qJ(2));
t84 = t78 * pkin(2) + pkin(5);
t80 = cos(qJ(2));
t67 = t80 * pkin(2) + pkin(1);
t77 = -pkin(6) - qJ(3);
t76 = qJ(2) + pkin(9);
t75 = -pkin(7) + t77;
t68 = sin(t76);
t83 = pkin(3) * t68 + t84;
t69 = cos(t76);
t58 = pkin(3) * t69 + t67;
t70 = qJ(4) + t76;
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t82 = g(1) * t81 + g(2) * t79;
t71 = -pkin(8) + t75;
t66 = qJ(5) + t70;
t65 = cos(t70);
t64 = sin(t70);
t61 = cos(t66);
t60 = sin(t66);
t59 = g(1) * t79 - g(2) * t81;
t57 = pkin(4) * t65 + t58;
t1 = [0, 0, 0, 0, 0, 0, -t82, t59, -g(3), -t85, 0, 0, 0, 0, 0, 0, -g(3) * t78 - t80 * t82, -g(3) * t80 + t78 * t82, -t59, -g(1) * (pkin(1) * t81 + pkin(6) * t79) - g(2) * (pkin(1) * t79 - pkin(6) * t81) - t85, 0, 0, 0, 0, 0, 0, -g(3) * t68 - t69 * t82, -g(3) * t69 + t68 * t82, -t59, -g(1) * (t67 * t81 - t77 * t79) - g(2) * (t67 * t79 + t77 * t81) - g(3) * t84, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t65 * t82, -g(3) * t65 + t64 * t82, -t59, -g(1) * (t58 * t81 - t75 * t79) - g(2) * (t58 * t79 + t75 * t81) - g(3) * t83, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t61 * t82, -g(3) * t61 + t60 * t82, -t59, -g(1) * (t57 * t81 - t71 * t79) - g(2) * (t57 * t79 + t71 * t81) - g(3) * (pkin(4) * t64 + t83);];
U_reg = t1;

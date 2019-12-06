% Calculate inertial parameters regressor of potential energy for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:50
% EndTime: 2019-12-05 18:38:50
% DurationCPUTime: 0.09s
% Computational Cost: add. (118->45), mult. (99->54), div. (0->0), fcn. (86->10), ass. (0->25)
t85 = g(3) * pkin(5);
t81 = -pkin(7) - pkin(6);
t77 = sin(qJ(2));
t84 = t77 * pkin(2) + pkin(5);
t79 = cos(qJ(2));
t67 = t79 * pkin(2) + pkin(1);
t76 = qJ(2) + qJ(3);
t75 = -qJ(4) + t81;
t69 = sin(t76);
t83 = pkin(3) * t69 + t84;
t70 = cos(t76);
t58 = pkin(3) * t70 + t67;
t68 = pkin(9) + t76;
t78 = sin(qJ(1));
t80 = cos(qJ(1));
t82 = g(1) * t80 + g(2) * t78;
t71 = -pkin(8) + t75;
t66 = qJ(5) + t68;
t63 = cos(t68);
t62 = sin(t68);
t61 = cos(t66);
t60 = sin(t66);
t59 = g(1) * t78 - g(2) * t80;
t57 = pkin(4) * t63 + t58;
t1 = [0, 0, 0, 0, 0, 0, -t82, t59, -g(3), -t85, 0, 0, 0, 0, 0, 0, -g(3) * t77 - t82 * t79, -g(3) * t79 + t82 * t77, -t59, -g(1) * (t80 * pkin(1) + t78 * pkin(6)) - g(2) * (t78 * pkin(1) - t80 * pkin(6)) - t85, 0, 0, 0, 0, 0, 0, -g(3) * t69 - t82 * t70, -g(3) * t70 + t82 * t69, -t59, -g(1) * (t80 * t67 - t78 * t81) - g(2) * (t78 * t67 + t80 * t81) - g(3) * t84, 0, 0, 0, 0, 0, 0, -g(3) * t62 - t82 * t63, -g(3) * t63 + t82 * t62, -t59, -g(1) * (t80 * t58 - t78 * t75) - g(2) * (t78 * t58 + t80 * t75) - g(3) * t83, 0, 0, 0, 0, 0, 0, -g(3) * t60 - t82 * t61, -g(3) * t61 + t82 * t60, -t59, -g(1) * (t80 * t57 - t78 * t71) - g(2) * (t78 * t57 + t80 * t71) - g(3) * (pkin(4) * t62 + t83);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:16
% EndTime: 2019-12-05 16:49:17
% DurationCPUTime: 0.17s
% Computational Cost: add. (122->66), mult. (169->86), div. (0->0), fcn. (168->8), ass. (0->34)
t66 = -pkin(7) - pkin(6);
t63 = sin(qJ(2));
t80 = g(3) * t63;
t62 = sin(qJ(3));
t79 = t62 * pkin(3);
t64 = cos(qJ(3));
t50 = t64 * pkin(3) + pkin(2);
t78 = g(3) * qJ(1);
t58 = -qJ(5) + t66;
t77 = t58 * t63;
t60 = sin(pkin(8));
t76 = t60 * t62;
t65 = cos(qJ(2));
t75 = t60 * t65;
t61 = cos(pkin(8));
t74 = t61 * t65;
t73 = t62 * t65;
t72 = t63 * t66;
t71 = t64 * t65;
t70 = t61 * pkin(1) + t60 * pkin(5);
t54 = t60 * pkin(1);
t69 = -t61 * pkin(5) + t54;
t68 = pkin(2) * t65 + pkin(6) * t63;
t67 = g(1) * t61 + g(2) * t60;
t59 = qJ(3) + qJ(4);
t52 = cos(t59);
t51 = sin(t59);
t49 = g(1) * t60 - g(2) * t61;
t48 = pkin(4) * t51 + t79;
t47 = pkin(4) * t52 + t50;
t46 = -g(3) * t65 + t63 * t67;
t45 = -g(1) * (t60 * t51 + t52 * t74) - g(2) * (-t61 * t51 + t52 * t75) - t52 * t80;
t44 = -g(1) * (-t51 * t74 + t60 * t52) - g(2) * (-t51 * t75 - t61 * t52) + t51 * t80;
t1 = [0, 0, 0, 0, 0, 0, -t67, t49, -g(3), -t78, 0, 0, 0, 0, 0, 0, -t65 * t67 - t80, t46, -t49, -g(1) * t70 - g(2) * t69 - t78, 0, 0, 0, 0, 0, 0, -g(1) * (t61 * t71 + t76) - g(2) * (t60 * t71 - t61 * t62) - t64 * t80, -g(1) * (t60 * t64 - t61 * t73) - g(2) * (-t60 * t73 - t61 * t64) + t62 * t80, -t46, -g(1) * (t61 * t68 + t70) - g(2) * (t60 * t68 + t69) - g(3) * (t63 * pkin(2) - t65 * pkin(6) + qJ(1)), 0, 0, 0, 0, 0, 0, t45, t44, -t46, -g(1) * (pkin(3) * t76 + t70) - g(2) * (t50 * t75 - t60 * t72 + t54) - g(3) * (t63 * t50 + t65 * t66 + qJ(1)) + (-g(1) * (t50 * t65 - t72) - g(2) * (-pkin(5) - t79)) * t61, 0, 0, 0, 0, 0, 0, t45, t44, -t46, -g(1) * (t60 * t48 + t70) - g(2) * (t47 * t75 - t60 * t77 + t54) - g(3) * (t63 * t47 + t65 * t58 + qJ(1)) + (-g(1) * (t47 * t65 - t77) - g(2) * (-pkin(5) - t48)) * t61;];
U_reg = t1;

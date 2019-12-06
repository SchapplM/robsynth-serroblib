% Calculate inertial parameters regressor of potential energy for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:15
% EndTime: 2019-12-05 15:11:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (125->57), mult. (268->78), div. (0->0), fcn. (301->8), ass. (0->43)
t65 = sin(pkin(8));
t67 = cos(pkin(8));
t95 = pkin(2) * t67 + pkin(5) * t65;
t92 = g(3) * qJ(1);
t69 = sin(qJ(4));
t91 = t65 * t69;
t70 = sin(qJ(3));
t90 = t65 * t70;
t71 = cos(qJ(4));
t89 = t65 * t71;
t72 = cos(qJ(3));
t88 = t65 * t72;
t66 = sin(pkin(7));
t87 = t66 * t70;
t86 = t66 * t72;
t68 = cos(pkin(7));
t85 = t68 * t70;
t84 = t68 * t72;
t83 = t68 * pkin(1) + t66 * qJ(2);
t82 = t66 * pkin(1) - t68 * qJ(2);
t81 = t95 * t68 + t83;
t80 = t65 * pkin(2) - t67 * pkin(5) + qJ(1);
t79 = g(1) * t68 + g(2) * t66;
t78 = t95 * t66 + t82;
t77 = pkin(3) * t88 + pkin(6) * t90 + t80;
t46 = t67 * t86 - t85;
t38 = t46 * t69 - t66 * t89;
t48 = t67 * t84 + t87;
t40 = t48 * t69 - t68 * t89;
t49 = t67 * t71 + t69 * t88;
t76 = g(1) * t40 + g(2) * t38 + g(3) * t49;
t47 = t67 * t85 - t86;
t75 = t48 * pkin(3) + t47 * pkin(6) + t81;
t45 = t67 * t87 + t84;
t74 = g(1) * t47 + g(2) * t45 + g(3) * t90;
t73 = t46 * pkin(3) + t45 * pkin(6) + t78;
t51 = g(1) * t66 - g(2) * t68;
t50 = -t67 * t69 + t71 * t88;
t42 = -g(3) * t67 + t79 * t65;
t41 = t48 * t71 + t68 * t91;
t39 = t46 * t71 + t66 * t91;
t36 = -g(1) * t41 - g(2) * t39 - g(3) * t50;
t1 = [0, 0, 0, 0, 0, 0, -t79, t51, -g(3), -t92, 0, 0, 0, 0, 0, 0, -g(3) * t65 - t79 * t67, t42, -t51, -g(1) * t83 - g(2) * t82 - t92, 0, 0, 0, 0, 0, 0, -g(1) * t48 - g(2) * t46 - g(3) * t88, t74, -t42, -g(1) * t81 - g(2) * t78 - g(3) * t80, 0, 0, 0, 0, 0, 0, t36, t76, -t74, -g(1) * t75 - g(2) * t73 - g(3) * t77, 0, 0, 0, 0, 0, 0, t36, -t74, -t76, -g(1) * (t41 * pkin(4) + t40 * qJ(5) + t75) - g(2) * (t39 * pkin(4) + t38 * qJ(5) + t73) - g(3) * (t50 * pkin(4) + t49 * qJ(5) + t77);];
U_reg = t1;

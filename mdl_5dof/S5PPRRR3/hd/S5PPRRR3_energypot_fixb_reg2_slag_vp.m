% Calculate inertial parameters regressor of potential energy for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:58
% EndTime: 2019-12-05 15:16:58
% DurationCPUTime: 0.20s
% Computational Cost: add. (130->74), mult. (243->105), div. (0->0), fcn. (266->10), ass. (0->42)
t62 = cos(pkin(9));
t90 = pkin(2) * t62;
t60 = sin(pkin(9));
t89 = g(3) * t60;
t88 = g(3) * qJ(1);
t61 = sin(pkin(8));
t87 = t60 * t61;
t63 = cos(pkin(8));
t86 = t60 * t63;
t64 = sin(qJ(4));
t85 = t60 * t64;
t66 = cos(qJ(4));
t84 = t60 * t66;
t67 = cos(qJ(3));
t83 = t60 * t67;
t65 = sin(qJ(3));
t82 = t61 * t65;
t81 = t61 * t67;
t80 = t63 * t65;
t79 = t63 * t67;
t78 = t63 * pkin(1) + t61 * qJ(2);
t77 = t60 * pkin(2) + qJ(1);
t76 = t61 * t85;
t75 = t63 * t85;
t74 = t61 * pkin(1) - t63 * qJ(2);
t73 = pkin(5) * t86 + t63 * t90 + t78;
t72 = -t62 * pkin(5) + t77;
t71 = g(1) * t63 + g(2) * t61;
t70 = pkin(5) * t87 + t61 * t90 + t74;
t42 = t62 * t82 + t79;
t44 = t62 * t80 - t81;
t69 = g(1) * t44 + g(2) * t42 + t65 * t89;
t68 = -pkin(7) - pkin(6);
t59 = qJ(4) + qJ(5);
t55 = cos(t59);
t54 = sin(t59);
t52 = t66 * pkin(4) + pkin(3);
t46 = g(1) * t61 - g(2) * t63;
t45 = t62 * t79 + t82;
t43 = t62 * t81 - t80;
t41 = -g(3) * t62 + t71 * t60;
t1 = [0, 0, 0, 0, 0, 0, -t71, t46, -g(3), -t88, 0, 0, 0, 0, 0, 0, -t62 * t71 - t89, t41, -t46, -g(1) * t78 - g(2) * t74 - t88, 0, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t43 - g(3) * t83, t69, -t41, -g(1) * t73 - g(2) * t70 - g(3) * t72, 0, 0, 0, 0, 0, 0, -g(1) * (t45 * t66 + t75) - g(2) * (t43 * t66 + t76) - g(3) * (-t62 * t64 + t66 * t83), -g(1) * (-t45 * t64 + t63 * t84) - g(2) * (-t43 * t64 + t61 * t84) - g(3) * (-t62 * t66 - t64 * t83), -t69, -g(1) * (t45 * pkin(3) + t44 * pkin(6) + t73) - g(2) * (t43 * pkin(3) + t42 * pkin(6) + t70) - g(3) * ((pkin(3) * t67 + pkin(6) * t65) * t60 + t72), 0, 0, 0, 0, 0, 0, -g(1) * (t45 * t55 + t54 * t86) - g(2) * (t43 * t55 + t54 * t87) - g(3) * (-t62 * t54 + t55 * t83), -g(1) * (-t45 * t54 + t55 * t86) - g(2) * (-t43 * t54 + t55 * t87) - g(3) * (-t54 * t83 - t62 * t55), -t69, -g(1) * (pkin(4) * t75 - t44 * t68 + t45 * t52 + t73) - g(2) * (pkin(4) * t76 - t42 * t68 + t43 * t52 + t70) - g(3) * ((-pkin(4) * t64 - pkin(5)) * t62 + (t52 * t67 - t65 * t68) * t60 + t77);];
U_reg = t1;

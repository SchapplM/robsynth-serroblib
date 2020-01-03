% Calculate inertial parameters regressor of potential energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:47
% EndTime: 2019-12-31 16:36:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (114->56), mult. (262->90), div. (0->0), fcn. (311->10), ass. (0->37)
t58 = sin(pkin(4));
t82 = pkin(5) * t58;
t62 = sin(qJ(3));
t81 = t58 * t62;
t63 = sin(qJ(2));
t80 = t58 * t63;
t65 = cos(qJ(3));
t79 = t58 * t65;
t66 = cos(qJ(2));
t78 = t58 * t66;
t60 = cos(pkin(4));
t77 = t60 * t63;
t76 = t60 * t66;
t57 = sin(pkin(8));
t59 = cos(pkin(8));
t75 = t59 * pkin(1) + t57 * t82;
t74 = t60 * pkin(5) + qJ(1);
t73 = t57 * pkin(1) - t59 * t82;
t72 = g(1) * t57 - g(2) * t59;
t45 = t57 * t76 + t59 * t63;
t46 = -t57 * t77 + t59 * t66;
t71 = t46 * pkin(2) + t45 * pkin(6) + t75;
t70 = pkin(2) * t80 - pkin(6) * t78 + t74;
t44 = t57 * t66 + t59 * t77;
t37 = t44 * t62 + t59 * t79;
t39 = t46 * t62 - t57 * t79;
t47 = -t60 * t65 + t62 * t80;
t69 = g(1) * t39 + g(2) * t37 + g(3) * t47;
t43 = t57 * t63 - t59 * t76;
t68 = t44 * pkin(2) + t43 * pkin(6) + t73;
t67 = -g(1) * t45 - g(2) * t43 + g(3) * t78;
t64 = cos(qJ(4));
t61 = sin(qJ(4));
t48 = t60 * t62 + t63 * t79;
t40 = t46 * t65 + t57 * t81;
t38 = t44 * t65 - t59 * t81;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t59 - g(2) * t57, t72, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t46 - g(2) * t44 - g(3) * t80, -t67, -g(3) * t60 - t72 * t58, -g(1) * t75 - g(2) * t73 - g(3) * t74, 0, 0, 0, 0, 0, 0, -g(1) * t40 - g(2) * t38 - g(3) * t48, t69, t67, -g(1) * t71 - g(2) * t68 - g(3) * t70, 0, 0, 0, 0, 0, 0, -g(1) * (t40 * t64 + t45 * t61) - g(2) * (t38 * t64 + t43 * t61) - g(3) * (t48 * t64 - t61 * t78), -g(1) * (-t40 * t61 + t45 * t64) - g(2) * (-t38 * t61 + t43 * t64) - g(3) * (-t48 * t61 - t64 * t78), -t69, -g(1) * (t40 * pkin(3) + t39 * pkin(7) + t71) - g(2) * (t38 * pkin(3) + t37 * pkin(7) + t68) - g(3) * (t48 * pkin(3) + t47 * pkin(7) + t70);];
U_reg = t1;

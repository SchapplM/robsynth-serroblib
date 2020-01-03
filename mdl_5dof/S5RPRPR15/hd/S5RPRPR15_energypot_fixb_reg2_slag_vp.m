% Calculate inertial parameters regressor of potential energy for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR15_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:23
% EndTime: 2019-12-31 18:37:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (87->58), mult. (132->71), div. (0->0), fcn. (127->8), ass. (0->35)
t56 = cos(pkin(8));
t46 = t56 * pkin(4) + pkin(3);
t57 = -pkin(7) - qJ(4);
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t81 = t46 * t58 + t57 * t60;
t80 = g(3) * pkin(5);
t79 = pkin(2) + pkin(5);
t61 = cos(qJ(1));
t78 = g(2) * t61;
t77 = g(3) * t60;
t54 = pkin(8) + qJ(5);
t47 = sin(t54);
t59 = sin(qJ(1));
t74 = t59 * t47;
t48 = cos(t54);
t73 = t59 * t48;
t55 = sin(pkin(8));
t72 = t59 * t55;
t71 = t59 * t56;
t70 = t61 * t47;
t69 = t61 * t48;
t68 = t61 * t55;
t67 = t61 * t56;
t50 = t59 * pkin(6);
t51 = t59 * pkin(1);
t66 = t50 + t51;
t65 = t61 * pkin(1) + t59 * qJ(2);
t64 = t61 * pkin(6) + t65;
t63 = -t61 * qJ(2) + t51;
t44 = g(1) * t59 - t78;
t62 = pkin(3) * t58 - qJ(4) * t60;
t45 = g(1) * t61 + g(2) * t59;
t43 = -g(3) * t58 + t44 * t60;
t1 = [0, 0, 0, 0, 0, 0, -t45, t44, -g(3), -t80, 0, 0, 0, 0, 0, 0, -g(3), t45, -t44, -g(1) * t65 - g(2) * t63 - t80, 0, 0, 0, 0, 0, 0, -t44 * t58 - t77, -t43, -t45, -g(1) * t64 - g(2) * (t50 + t63) - g(3) * t79, 0, 0, 0, 0, 0, 0, -g(1) * (t58 * t71 + t68) - g(2) * (-t58 * t67 + t72) - t56 * t77, -g(1) * (-t58 * t72 + t67) - g(2) * (t58 * t68 + t71) + t55 * t77, t43, -g(1) * (t62 * t59 + t64) - g(2) * t66 - g(3) * (t60 * pkin(3) + t58 * qJ(4) + t79) - (-qJ(2) - t62) * t78, 0, 0, 0, 0, 0, 0, -g(1) * (t58 * t73 + t70) - g(2) * (-t58 * t69 + t74) - t48 * t77, -g(1) * (-t58 * t74 + t69) - g(2) * (t58 * t70 + t73) + t47 * t77, t43, -g(1) * (t81 * t59 + t64) - g(2) * (pkin(4) * t72 + t66) - g(3) * (t60 * t46 - t58 * t57 + t79) + (-g(1) * pkin(4) * t55 - g(2) * (-qJ(2) - t81)) * t61;];
U_reg = t1;

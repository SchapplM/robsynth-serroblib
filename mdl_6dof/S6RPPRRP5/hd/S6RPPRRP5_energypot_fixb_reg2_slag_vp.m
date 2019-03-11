% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:46
% EndTime: 2019-03-09 02:08:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (96->57), mult. (151->62), div. (0->0), fcn. (143->6), ass. (0->35)
t82 = g(3) * pkin(6);
t81 = pkin(2) + pkin(6);
t59 = sin(qJ(5));
t80 = pkin(5) * t59;
t63 = cos(qJ(4));
t79 = g(3) * t63;
t61 = sin(qJ(1));
t78 = t61 * t59;
t62 = cos(qJ(5));
t77 = t61 * t62;
t64 = cos(qJ(1));
t76 = t64 * t59;
t75 = t64 * t62;
t74 = t64 * pkin(1) + t61 * qJ(2);
t73 = pkin(3) + t81;
t72 = t64 * qJ(3) + t74;
t54 = t61 * pkin(1);
t71 = -t64 * qJ(2) + t54;
t70 = g(1) * t72;
t51 = t61 * qJ(3);
t69 = t51 + t71;
t60 = sin(qJ(4));
t68 = pkin(4) * t60 - pkin(8) * t63;
t49 = g(1) * t64 + g(2) * t61;
t50 = t62 * pkin(5) + pkin(4);
t58 = -qJ(6) - pkin(8);
t67 = t50 * t60 + t58 * t63;
t66 = -t61 * pkin(7) + t72;
t55 = t64 * pkin(7);
t65 = t55 + t69;
t48 = g(1) * t61 - g(2) * t64;
t47 = -g(3) * t60 + t49 * t63;
t46 = -g(1) * (t60 * t75 - t78) - g(2) * (t60 * t77 + t76) - t62 * t79;
t45 = -g(1) * (-t60 * t76 - t77) - g(2) * (-t60 * t78 + t75) + t59 * t79;
t1 = [0, 0, 0, 0, 0, 0, -t49, t48, -g(3), -t82, 0, 0, 0, 0, 0, 0, -g(3), t49, -t48, -g(1) * t74 - g(2) * t71 - t82, 0, 0, 0, 0, 0, 0, -g(3), -t48, -t49, -g(2) * t69 - g(3) * t81 - t70, 0, 0, 0, 0, 0, 0, -t49 * t60 - t79, -t47, t48, -g(1) * t66 - g(2) * t65 - g(3) * t73, 0, 0, 0, 0, 0, 0, t46, t45, t47, -g(1) * (t68 * t64 + t66) - g(2) * (t68 * t61 + t65) - g(3) * (t63 * pkin(4) + t60 * pkin(8) + t73) 0, 0, 0, 0, 0, 0, t46, t45, t47, -t70 - g(2) * (t51 + t54 + t55) - g(3) * (t63 * t50 - t60 * t58 + t73) + (-g(1) * t67 - g(2) * (-qJ(2) + t80)) * t64 + (-g(1) * (-pkin(7) - t80) - g(2) * t67) * t61;];
U_reg  = t1;

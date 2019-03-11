% Calculate inertial parameters regressor of potential energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:01
% EndTime: 2019-03-09 01:32:01
% DurationCPUTime: 0.11s
% Computational Cost: add. (167->58), mult. (129->66), div. (0->0), fcn. (117->10), ass. (0->35)
t64 = pkin(10) + qJ(5);
t57 = sin(t64);
t59 = cos(t64);
t92 = pkin(5) * t57 - pkin(8) * t59;
t66 = sin(pkin(10));
t91 = pkin(4) * t66;
t88 = g(3) * t59;
t68 = qJ(2) + pkin(6);
t87 = g(3) * t68;
t65 = qJ(1) + pkin(9);
t58 = sin(t65);
t70 = sin(qJ(6));
t86 = t58 * t70;
t72 = cos(qJ(6));
t85 = t58 * t72;
t60 = cos(t65);
t84 = t60 * t70;
t83 = t60 * t72;
t71 = sin(qJ(1));
t82 = t71 * pkin(1) + t58 * pkin(2);
t81 = pkin(3) + t68;
t73 = cos(qJ(1));
t80 = t73 * pkin(1) + t60 * pkin(2) + t58 * qJ(3);
t79 = -qJ(3) - t91;
t67 = cos(pkin(10));
t78 = t67 * pkin(4) + t81;
t77 = t58 * t91 + t80;
t69 = -pkin(7) - qJ(4);
t76 = -t58 * t69 + t82;
t75 = -t60 * qJ(3) + t82;
t51 = g(1) * t58 - g(2) * t60;
t74 = -g(1) * t73 - g(2) * t71;
t52 = g(1) * t60 + g(2) * t58;
t50 = -g(3) * t57 + t51 * t59;
t1 = [0, 0, 0, 0, 0, 0, t74, g(1) * t71 - g(2) * t73, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t52, t51, -g(3), t74 * pkin(1) - t87, 0, 0, 0, 0, 0, 0, -g(3), t52, -t51, -g(1) * t80 - g(2) * t75 - t87, 0, 0, 0, 0, 0, 0, -g(3) * t67 - t51 * t66, g(3) * t66 - t51 * t67, -t52, -g(1) * (t60 * qJ(4) + t80) - g(2) * (t58 * qJ(4) + t75) - g(3) * t81, 0, 0, 0, 0, 0, 0, -t51 * t57 - t88, -t50, -t52, -g(1) * (-t60 * t69 + t77) - g(2) * (t79 * t60 + t76) - g(3) * t78, 0, 0, 0, 0, 0, 0, -g(1) * (t57 * t85 + t84) - g(2) * (-t57 * t83 + t86) - t72 * t88, -g(1) * (-t57 * t86 + t83) - g(2) * (t57 * t84 + t85) + t70 * t88, t50, -g(1) * (t92 * t58 + t77) - g(2) * t76 - g(3) * (t59 * pkin(5) + t57 * pkin(8) + t78) + (g(1) * t69 - g(2) * (t79 - t92)) * t60;];
U_reg  = t1;

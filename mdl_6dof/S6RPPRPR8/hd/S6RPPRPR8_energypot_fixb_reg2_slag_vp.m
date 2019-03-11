% Calculate inertial parameters regressor of potential energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:09
% EndTime: 2019-03-09 01:56:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (132->66), mult. (155->73), div. (0->0), fcn. (143->8), ass. (0->39)
t100 = g(3) * pkin(6);
t99 = pkin(2) + pkin(6);
t72 = sin(pkin(9));
t98 = pkin(3) * t72;
t71 = pkin(9) + qJ(4);
t65 = sin(t71);
t97 = pkin(8) * t65;
t78 = cos(qJ(1));
t96 = g(2) * t78;
t95 = g(3) * t65;
t74 = -pkin(7) - qJ(3);
t94 = pkin(5) - t74;
t76 = sin(qJ(1));
t93 = t65 * t76;
t75 = sin(qJ(6));
t92 = t76 * t75;
t77 = cos(qJ(6));
t91 = t76 * t77;
t90 = t78 * t75;
t89 = t78 * t77;
t88 = t78 * pkin(1) + t76 * qJ(2);
t66 = cos(t71);
t87 = qJ(5) * t66;
t73 = cos(pkin(9));
t86 = t73 * pkin(3) + t99;
t85 = t76 * t98 + t88;
t84 = -qJ(2) - t98;
t69 = t76 * pkin(1);
t83 = -t76 * t74 + t69;
t82 = -t78 * qJ(2) + t69;
t81 = pkin(4) * t93 + t85;
t80 = t66 * pkin(4) + t65 * qJ(5) + t86;
t58 = g(1) * t76 - t96;
t79 = -pkin(4) * t65 + t84;
t59 = g(1) * t78 + g(2) * t76;
t57 = t78 * t87;
t56 = t58 * t66 - t95;
t55 = g(1) * t93 + g(3) * t66 - t65 * t96;
t1 = [0, 0, 0, 0, 0, 0, -t59, t58, -g(3), -t100, 0, 0, 0, 0, 0, 0, -g(3), t59, -t58, -g(1) * t88 - g(2) * t82 - t100, 0, 0, 0, 0, 0, 0, -g(3) * t73 - t58 * t72, g(3) * t72 - t58 * t73, -t59, -g(1) * (t78 * qJ(3) + t88) - g(2) * (t76 * qJ(3) + t82) - g(3) * t99, 0, 0, 0, 0, 0, 0, -t55, -t56, -t59, -g(1) * (-t78 * t74 + t85) - g(2) * (t84 * t78 + t83) - g(3) * t86, 0, 0, 0, 0, 0, 0, -t59, t55, t56, -g(1) * (-t76 * t87 + t81) - g(2) * (t57 + t83) - g(3) * t80 + (g(1) * t74 - g(2) * t79) * t78, 0, 0, 0, 0, 0, 0, -g(1) * (-t66 * t92 + t89) - g(2) * (t66 * t90 + t91) - t75 * t95, -g(1) * (-t66 * t91 - t90) - g(2) * (t66 * t89 - t92) - t77 * t95, -t55, -g(1) * t81 - g(2) * (t57 + t69) - g(3) * (t66 * pkin(8) + t80) + (-g(1) * (-t87 + t97) - g(2) * t94) * t76 + (-g(1) * t94 - g(2) * (t79 - t97)) * t78;];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:58
% EndTime: 2019-03-09 01:33:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (152->57), mult. (217->72), div. (0->0), fcn. (247->10), ass. (0->35)
t98 = g(3) * pkin(6);
t73 = pkin(10) + qJ(5);
t66 = sin(t73);
t97 = g(3) * t66;
t76 = -qJ(3) + pkin(6);
t96 = g(3) * t76;
t95 = cos(qJ(1));
t94 = sin(qJ(1));
t67 = cos(t73);
t78 = sin(qJ(6));
t93 = t67 * t78;
t79 = cos(qJ(6));
t92 = t67 * t79;
t91 = t95 * pkin(1) + t94 * qJ(2);
t90 = cos(pkin(9));
t89 = sin(pkin(9));
t88 = t95 * pkin(2) + t91;
t74 = sin(pkin(10));
t87 = -t74 * pkin(4) + t76;
t58 = -t94 * t89 - t95 * t90;
t59 = t95 * t89 - t94 * t90;
t75 = cos(pkin(10));
t65 = t75 * pkin(4) + pkin(3);
t77 = -pkin(7) - qJ(4);
t86 = -t58 * t65 - t59 * t77 + t88;
t85 = -pkin(5) * t67 - pkin(8) * t66;
t84 = g(1) * t59 - g(2) * t58;
t83 = g(1) * t58 + g(2) * t59;
t82 = t94 * pkin(1) - t95 * qJ(2);
t81 = t94 * pkin(2) + t82;
t80 = t58 * t77 - t59 * t65 + t81;
t61 = -g(1) * t95 - g(2) * t94;
t60 = g(1) * t94 - g(2) * t95;
t52 = g(3) * t67 - t83 * t66;
t1 = [0, 0, 0, 0, 0, 0, t61, t60, -g(3), -t98, 0, 0, 0, 0, 0, 0, t61, -g(3), -t60, -g(1) * t91 - g(2) * t82 - t98, 0, 0, 0, 0, 0, 0, t83, t84, g(3), -g(1) * t88 - g(2) * t81 - t96, 0, 0, 0, 0, 0, 0, g(3) * t74 + t83 * t75, g(3) * t75 - t83 * t74, -t84, -g(1) * (-t58 * pkin(3) + t59 * qJ(4) + t88) - g(2) * (-t59 * pkin(3) - t58 * qJ(4) + t81) - t96, 0, 0, 0, 0, 0, 0, t83 * t67 + t97, t52, -t84, -g(1) * t86 - g(2) * t80 - g(3) * t87, 0, 0, 0, 0, 0, 0, -g(1) * (-t58 * t92 + t59 * t78) - g(2) * (-t58 * t78 - t59 * t92) + t79 * t97, -g(1) * (t58 * t93 + t59 * t79) - g(2) * (-t58 * t79 + t59 * t93) - t78 * t97, -t52, -g(1) * (t85 * t58 + t86) - g(2) * (t85 * t59 + t80) - g(3) * (-t66 * pkin(5) + t67 * pkin(8) + t87);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:02
% EndTime: 2019-03-09 03:03:02
% DurationCPUTime: 0.14s
% Computational Cost: add. (201->62), mult. (160->73), div. (0->0), fcn. (152->10), ass. (0->38)
t86 = cos(qJ(5));
t69 = t86 * pkin(5) + pkin(4);
t78 = qJ(3) + pkin(10);
t71 = sin(t78);
t73 = cos(t78);
t80 = -qJ(6) - pkin(8);
t104 = t69 * t73 - t71 * t80;
t103 = g(3) * t71;
t82 = qJ(2) + pkin(6);
t102 = g(3) * t82;
t79 = qJ(1) + pkin(9);
t72 = sin(t79);
t83 = sin(qJ(5));
t99 = t72 * t83;
t98 = t72 * t86;
t74 = cos(t79);
t97 = t74 * t83;
t96 = t74 * t86;
t87 = cos(qJ(3));
t70 = t87 * pkin(3) + pkin(2);
t88 = cos(qJ(1));
t77 = t88 * pkin(1);
t95 = t74 * t70 + t77;
t85 = sin(qJ(1));
t76 = t85 * pkin(1);
t81 = -qJ(4) - pkin(7);
t94 = t72 * t70 + t74 * t81 + t76;
t84 = sin(qJ(3));
t93 = t84 * pkin(3) + t82;
t92 = -t72 * t81 + t95;
t91 = pkin(4) * t73 + pkin(8) * t71;
t90 = g(1) * t74 + g(2) * t72;
t89 = -g(1) * t88 - g(2) * t85;
t63 = g(1) * t72 - g(2) * t74;
t62 = -g(3) * t73 + t90 * t71;
t61 = -g(1) * (t73 * t96 + t99) - g(2) * (t73 * t98 - t97) - t86 * t103;
t60 = -g(1) * (-t73 * t97 + t98) - g(2) * (-t73 * t99 - t96) + t83 * t103;
t1 = [0, 0, 0, 0, 0, 0, t89, g(1) * t85 - g(2) * t88, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t90, t63, -g(3), t89 * pkin(1) - t102, 0, 0, 0, 0, 0, 0, -g(3) * t84 - t90 * t87, -g(3) * t87 + t90 * t84, -t63, -g(1) * (t74 * pkin(2) + t72 * pkin(7) + t77) - g(2) * (t72 * pkin(2) - t74 * pkin(7) + t76) - t102, 0, 0, 0, 0, 0, 0, -t90 * t73 - t103, t62, -t63, -g(1) * t92 - g(2) * t94 - g(3) * t93, 0, 0, 0, 0, 0, 0, t61, t60, -t62, -g(1) * (t91 * t74 + t92) - g(2) * (t91 * t72 + t94) - g(3) * (t71 * pkin(4) - t73 * pkin(8) + t93) 0, 0, 0, 0, 0, 0, t61, t60, -t62, -g(1) * (t104 * t74 + t95) - g(2) * (-pkin(5) * t97 + t94) - g(3) * (t71 * t69 + t73 * t80 + t93) + (-g(1) * (pkin(5) * t83 - t81) - g(2) * t104) * t72;];
U_reg  = t1;

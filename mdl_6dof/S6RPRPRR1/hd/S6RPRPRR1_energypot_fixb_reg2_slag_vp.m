% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:18
% EndTime: 2019-03-09 03:35:18
% DurationCPUTime: 0.11s
% Computational Cost: add. (200->62), mult. (136->72), div. (0->0), fcn. (124->12), ass. (0->38)
t83 = qJ(3) + pkin(11);
t77 = qJ(5) + t83;
t70 = sin(t77);
t105 = g(3) * t70;
t86 = qJ(2) + pkin(6);
t104 = g(3) * t86;
t91 = cos(qJ(3));
t72 = t91 * pkin(3) + pkin(2);
t84 = qJ(1) + pkin(10);
t74 = sin(t84);
t87 = sin(qJ(6));
t103 = t74 * t87;
t90 = cos(qJ(6));
t102 = t74 * t90;
t76 = cos(t84);
t101 = t76 * t87;
t100 = t76 * t90;
t85 = -qJ(4) - pkin(7);
t75 = cos(t83);
t66 = pkin(4) * t75 + t72;
t89 = sin(qJ(1));
t79 = t89 * pkin(1);
t82 = -pkin(8) + t85;
t99 = t74 * t66 + t76 * t82 + t79;
t88 = sin(qJ(3));
t98 = t88 * pkin(3) + t86;
t73 = sin(t83);
t97 = pkin(4) * t73 + t98;
t92 = cos(qJ(1));
t81 = t92 * pkin(1);
t96 = t76 * t66 - t74 * t82 + t81;
t71 = cos(t77);
t95 = pkin(5) * t71 + pkin(9) * t70;
t94 = g(1) * t76 + g(2) * t74;
t93 = -g(1) * t92 - g(2) * t89;
t65 = g(1) * t74 - g(2) * t76;
t62 = -g(3) * t71 + t94 * t70;
t1 = [0, 0, 0, 0, 0, 0, t93, g(1) * t89 - g(2) * t92, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t94, t65, -g(3), t93 * pkin(1) - t104, 0, 0, 0, 0, 0, 0, -g(3) * t88 - t94 * t91, -g(3) * t91 + t94 * t88, -t65, -g(1) * (t76 * pkin(2) + t74 * pkin(7) + t81) - g(2) * (t74 * pkin(2) - t76 * pkin(7) + t79) - t104, 0, 0, 0, 0, 0, 0, -g(3) * t73 - t94 * t75, -g(3) * t75 + t94 * t73, -t65, -g(1) * (t76 * t72 - t74 * t85 + t81) - g(2) * (t74 * t72 + t76 * t85 + t79) - g(3) * t98, 0, 0, 0, 0, 0, 0, -t94 * t71 - t105, t62, -t65, -g(1) * t96 - g(2) * t99 - g(3) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t100 + t103) - g(2) * (t71 * t102 - t101) - t90 * t105, -g(1) * (-t71 * t101 + t102) - g(2) * (-t71 * t103 - t100) + t87 * t105, -t62, -g(1) * (t95 * t76 + t96) - g(2) * (t95 * t74 + t99) - g(3) * (t70 * pkin(5) - t71 * pkin(9) + t97);];
U_reg  = t1;

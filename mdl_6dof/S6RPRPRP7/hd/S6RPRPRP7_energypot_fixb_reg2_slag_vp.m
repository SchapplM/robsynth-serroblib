% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:45
% EndTime: 2019-03-09 03:22:45
% DurationCPUTime: 0.16s
% Computational Cost: add. (140->63), mult. (168->73), div. (0->0), fcn. (160->8), ass. (0->36)
t69 = qJ(3) + pkin(9);
t63 = sin(t69);
t64 = cos(t69);
t96 = pkin(4) * t63 - pkin(8) * t64;
t95 = g(3) * pkin(6);
t94 = pkin(2) + pkin(6);
t73 = sin(qJ(3));
t93 = pkin(3) * t73;
t90 = g(3) * t64;
t72 = sin(qJ(5));
t74 = sin(qJ(1));
t89 = t74 * t72;
t75 = cos(qJ(5));
t88 = t74 * t75;
t77 = cos(qJ(1));
t87 = t77 * t72;
t86 = t77 * t75;
t85 = t77 * pkin(1) + t74 * qJ(2);
t76 = cos(qJ(3));
t84 = t76 * pkin(3) + t94;
t83 = t74 * t93 + t85;
t71 = -qJ(4) - pkin(7);
t82 = pkin(5) * t72 - t71;
t81 = -qJ(2) - t93;
t66 = t74 * pkin(1);
t80 = -t74 * t71 + t66;
t79 = -t77 * qJ(2) + t66;
t59 = g(1) * t74 - g(2) * t77;
t62 = t75 * pkin(5) + pkin(4);
t70 = -qJ(6) - pkin(8);
t78 = t62 * t63 + t64 * t70;
t60 = g(1) * t77 + g(2) * t74;
t58 = -g(3) * t63 + t59 * t64;
t57 = -g(1) * (t63 * t88 + t87) - g(2) * (-t63 * t86 + t89) - t75 * t90;
t56 = -g(1) * (-t63 * t89 + t86) - g(2) * (t63 * t87 + t88) + t72 * t90;
t1 = [0, 0, 0, 0, 0, 0, -t60, t59, -g(3), -t95, 0, 0, 0, 0, 0, 0, -g(3), t60, -t59, -g(1) * t85 - g(2) * t79 - t95, 0, 0, 0, 0, 0, 0, -g(3) * t76 - t59 * t73, g(3) * t73 - t59 * t76, -t60, -g(1) * (t77 * pkin(7) + t85) - g(2) * (t74 * pkin(7) + t79) - g(3) * t94, 0, 0, 0, 0, 0, 0, -t59 * t63 - t90, -t58, -t60, -g(1) * (-t77 * t71 + t83) - g(2) * (t81 * t77 + t80) - g(3) * t84, 0, 0, 0, 0, 0, 0, t57, t56, t58, -g(1) * (t74 * t96 + t83) - g(2) * t80 - g(3) * (t64 * pkin(4) + t63 * pkin(8) + t84) + (g(1) * t71 - g(2) * (t81 - t96)) * t77, 0, 0, 0, 0, 0, 0, t57, t56, t58, -g(1) * t83 - g(2) * t66 - g(3) * (t64 * t62 - t63 * t70 + t84) + (-g(1) * t78 - g(2) * t82) * t74 + (-g(1) * t82 - g(2) * (-t78 + t81)) * t77;];
U_reg  = t1;

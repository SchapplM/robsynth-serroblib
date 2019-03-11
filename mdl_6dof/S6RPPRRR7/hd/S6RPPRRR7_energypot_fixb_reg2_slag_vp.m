% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:47
% EndTime: 2019-03-09 02:33:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (149->64), mult. (144->71), div. (0->0), fcn. (132->10), ass. (0->36)
t70 = pkin(10) + qJ(4);
t64 = qJ(5) + t70;
t60 = sin(t64);
t61 = cos(t64);
t95 = pkin(5) * t60 - pkin(9) * t61;
t94 = g(3) * pkin(6);
t93 = pkin(2) + pkin(6);
t90 = g(3) * t61;
t71 = sin(pkin(10));
t89 = t71 * pkin(3);
t74 = sin(qJ(6));
t75 = sin(qJ(1));
t88 = t75 * t74;
t76 = cos(qJ(6));
t87 = t75 * t76;
t77 = cos(qJ(1));
t86 = t77 * t74;
t85 = t77 * t76;
t73 = -pkin(7) - qJ(3);
t84 = t77 * pkin(1) + t75 * qJ(2);
t62 = sin(t70);
t56 = pkin(4) * t62 + t89;
t83 = -qJ(2) - t56;
t72 = cos(pkin(10));
t82 = t72 * pkin(3) + t93;
t81 = t75 * t56 + t84;
t63 = cos(t70);
t80 = pkin(4) * t63 + t82;
t67 = t75 * pkin(1);
t69 = -pkin(8) + t73;
t79 = -t75 * t69 + t67;
t78 = -t77 * qJ(2) + t67;
t57 = g(1) * t75 - g(2) * t77;
t58 = g(1) * t77 + g(2) * t75;
t54 = -g(3) * t60 + t57 * t61;
t1 = [0, 0, 0, 0, 0, 0, -t58, t57, -g(3), -t94, 0, 0, 0, 0, 0, 0, -g(3), t58, -t57, -g(1) * t84 - g(2) * t78 - t94, 0, 0, 0, 0, 0, 0, -g(3) * t72 - t57 * t71, g(3) * t71 - t57 * t72, -t58, -g(1) * (t77 * qJ(3) + t84) - g(2) * (t75 * qJ(3) + t78) - g(3) * t93, 0, 0, 0, 0, 0, 0, -g(3) * t63 - t57 * t62, g(3) * t62 - t57 * t63, -t58, -g(1) * (-t77 * t73 + t75 * t89 + t84) - g(2) * (-t75 * t73 + t67 + (-qJ(2) - t89) * t77) - g(3) * t82, 0, 0, 0, 0, 0, 0, -t57 * t60 - t90, -t54, -t58, -g(1) * (-t77 * t69 + t81) - g(2) * (t77 * t83 + t79) - g(3) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (t60 * t87 + t86) - g(2) * (-t60 * t85 + t88) - t76 * t90, -g(1) * (-t60 * t88 + t85) - g(2) * (t60 * t86 + t87) + t74 * t90, t54, -g(1) * (t95 * t75 + t81) - g(2) * t79 - g(3) * (t61 * pkin(5) + t60 * pkin(9) + t80) + (g(1) * t69 - g(2) * (t83 - t95)) * t77;];
U_reg  = t1;

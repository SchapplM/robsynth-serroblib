% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:00
% EndTime: 2019-03-09 02:29:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (108->57), mult. (131->62), div. (0->0), fcn. (119->8), ass. (0->34)
t62 = qJ(4) + qJ(5);
t53 = sin(t62);
t54 = cos(t62);
t89 = pkin(5) * t53 - pkin(9) * t54;
t88 = g(3) * pkin(6);
t87 = pkin(2) + pkin(6);
t64 = sin(qJ(4));
t86 = pkin(4) * t64;
t83 = g(3) * t54;
t63 = sin(qJ(6));
t65 = sin(qJ(1));
t82 = t65 * t63;
t66 = cos(qJ(6));
t81 = t65 * t66;
t68 = cos(qJ(1));
t80 = t68 * t63;
t79 = t68 * t66;
t55 = t65 * qJ(3);
t58 = t65 * pkin(1);
t78 = t55 + t58;
t77 = t68 * pkin(1) + t65 * qJ(2);
t69 = -pkin(8) - pkin(7);
t76 = -qJ(2) - t69;
t75 = pkin(3) + t87;
t74 = t65 * t86 + t78;
t73 = t68 * qJ(3) + t77;
t67 = cos(qJ(4));
t72 = t67 * pkin(4) + t75;
t71 = -t68 * qJ(2) + t58;
t49 = g(1) * t68 + g(2) * t65;
t70 = g(1) * (t65 * t69 + t68 * t86 + t73);
t48 = g(1) * t65 - g(2) * t68;
t47 = -g(3) * t53 + t49 * t54;
t1 = [0, 0, 0, 0, 0, 0, -t49, t48, -g(3), -t88, 0, 0, 0, 0, 0, 0, -g(3), t49, -t48, -g(1) * t77 - g(2) * t71 - t88, 0, 0, 0, 0, 0, 0, -g(3), -t48, -t49, -g(1) * t73 - g(2) * (t55 + t71) - g(3) * t87, 0, 0, 0, 0, 0, 0, -g(3) * t67 - t49 * t64, g(3) * t64 - t49 * t67, t48, -g(1) * (-t65 * pkin(7) + t73) - g(2) * ((pkin(7) - qJ(2)) * t68 + t78) - g(3) * t75, 0, 0, 0, 0, 0, 0, -t49 * t53 - t83, -t47, t48, -t70 - g(2) * (t76 * t68 + t74) - g(3) * t72, 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t79 - t82) - g(2) * (t53 * t81 + t80) - t66 * t83, -g(1) * (-t53 * t80 - t81) - g(2) * (-t53 * t82 + t79) + t63 * t83, t47, -t70 - g(2) * (t89 * t65 + t74) - g(3) * (t54 * pkin(5) + t53 * pkin(9) + t72) + (-g(1) * t89 - g(2) * t76) * t68;];
U_reg  = t1;

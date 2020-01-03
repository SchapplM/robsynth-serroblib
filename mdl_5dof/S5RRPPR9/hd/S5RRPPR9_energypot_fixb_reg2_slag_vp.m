% Calculate inertial parameters regressor of potential energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:47
% EndTime: 2019-12-31 19:41:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (79->50), mult. (148->58), div. (0->0), fcn. (139->6), ass. (0->30)
t87 = g(3) * pkin(5);
t67 = cos(qJ(2));
t86 = g(3) * t67;
t64 = sin(qJ(2));
t85 = t64 * pkin(2) + pkin(5);
t63 = sin(qJ(5));
t65 = sin(qJ(1));
t84 = t65 * t63;
t66 = cos(qJ(5));
t83 = t65 * t66;
t82 = t65 * t67;
t68 = cos(qJ(1));
t81 = t67 * t68;
t80 = t68 * t63;
t79 = t68 * t66;
t78 = t68 * pkin(1) + t65 * pkin(6);
t77 = qJ(3) * t64;
t76 = t65 * pkin(1) - t68 * pkin(6);
t75 = pkin(2) * t81 + t68 * t77 + t78;
t74 = -t67 * qJ(3) + t85;
t73 = pkin(4) * t64 + pkin(7) * t67;
t72 = g(1) * t68 + g(2) * t65;
t71 = pkin(2) * t82 + t65 * t77 + t76;
t70 = pkin(3) * t82 + t68 * qJ(4) + t71;
t69 = pkin(3) * t81 - t65 * qJ(4) + t75;
t56 = t64 * pkin(3);
t48 = g(1) * t65 - g(2) * t68;
t47 = g(3) * t64 + t72 * t67;
t46 = t72 * t64 - t86;
t1 = [0, 0, 0, 0, 0, 0, -t72, t48, -g(3), -t87, 0, 0, 0, 0, 0, 0, -t47, t46, -t48, -g(1) * t78 - g(2) * t76 - t87, 0, 0, 0, 0, 0, 0, -t47, -t48, -t46, -g(1) * t75 - g(2) * t71 - g(3) * t74, 0, 0, 0, 0, 0, 0, -t46, t47, t48, -g(1) * t69 - g(2) * t70 - g(3) * (t56 + t74), 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t79 - t84) - g(2) * (t64 * t83 + t80) + t66 * t86, -g(1) * (-t64 * t80 - t83) - g(2) * (-t64 * t84 + t79) - t63 * t86, -t47, -g(1) * (t73 * t68 + t69) - g(2) * (t73 * t65 + t70) - g(3) * (t64 * pkin(7) + t56 + (-pkin(4) - qJ(3)) * t67 + t85);];
U_reg = t1;

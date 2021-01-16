% Calculate minimal parameter regressor of potential energy for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:23
% EndTime: 2021-01-16 00:10:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (78->31), mult. (96->41), div. (0->0), fcn. (101->8), ass. (0->23)
t67 = qJ(2) + qJ(3);
t65 = sin(t67);
t83 = g(3) * t65;
t69 = sin(qJ(4));
t71 = sin(qJ(1));
t82 = t71 * t69;
t72 = cos(qJ(4));
t81 = t71 * t72;
t74 = cos(qJ(1));
t80 = t74 * t69;
t79 = t74 * t72;
t78 = pkin(4) * t69 + pkin(6) + pkin(7);
t77 = g(1) * t74 + g(2) * t71;
t63 = t72 * pkin(4) + pkin(3);
t66 = cos(t67);
t68 = -qJ(5) - pkin(8);
t73 = cos(qJ(2));
t76 = t73 * pkin(2) + t63 * t66 - t65 * t68 + pkin(1);
t70 = sin(qJ(2));
t62 = -g(3) * t66 + t77 * t65;
t61 = -g(1) * (t66 * t79 + t82) - g(2) * (t66 * t81 - t80) - t72 * t83;
t60 = -g(1) * (-t66 * t80 + t81) - g(2) * (-t66 * t82 - t79) + t69 * t83;
t1 = [0, -t77, g(1) * t71 - g(2) * t74, 0, 0, 0, 0, 0, -g(3) * t70 - t77 * t73, -g(3) * t73 + t77 * t70, 0, 0, 0, 0, 0, -t77 * t66 - t83, t62, 0, 0, 0, 0, 0, t61, t60, t61, t60, -t62, -g(3) * (t70 * pkin(2) + t65 * t63 + t66 * t68 + pkin(5)) + (-g(1) * t76 + g(2) * t78) * t74 + (-g(1) * t78 - g(2) * t76) * t71;];
U_reg = t1;

% Calculate minimal parameter regressor of potential energy for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:25
% EndTime: 2021-01-15 12:06:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (68->29), mult. (63->43), div. (0->0), fcn. (61->10), ass. (0->23)
t71 = qJ(3) + pkin(9);
t67 = sin(t71);
t87 = g(3) * t67;
t72 = qJ(1) + pkin(8);
t68 = sin(t72);
t74 = sin(qJ(5));
t86 = t68 * t74;
t77 = cos(qJ(5));
t85 = t68 * t77;
t70 = cos(t72);
t84 = t70 * t74;
t83 = t70 * t77;
t82 = qJ(2) + pkin(5);
t81 = g(1) * t70 + g(2) * t68;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t80 = -g(1) * t79 - g(2) * t76;
t78 = cos(qJ(3));
t75 = sin(qJ(3));
t73 = -qJ(4) - pkin(6);
t69 = cos(t71);
t66 = t78 * pkin(3) + pkin(2);
t1 = [0, t80, g(1) * t76 - g(2) * t79, t80 * pkin(1) - g(3) * t82, 0, 0, 0, 0, 0, -g(3) * t75 - t81 * t78, -g(3) * t78 + t81 * t75, -t81 * t69 - t87, -g(3) * t69 + t81 * t67, -g(1) * t68 + g(2) * t70, -g(1) * (t79 * pkin(1) + t70 * t66 - t68 * t73) - g(2) * (t76 * pkin(1) + t68 * t66 + t70 * t73) - g(3) * (t75 * pkin(3) + t82), 0, 0, 0, 0, 0, -g(1) * (t69 * t83 + t86) - g(2) * (t69 * t85 - t84) - t77 * t87, -g(1) * (-t69 * t84 + t85) - g(2) * (-t69 * t86 - t83) + t74 * t87;];
U_reg = t1;

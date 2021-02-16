% Calculate minimal parameter regressor of potential energy for
% S5RPRPR5
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
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:55:11
% EndTime: 2021-01-15 11:55:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (83->46), mult. (104->70), div. (0->0), fcn. (110->10), ass. (0->32)
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t77 = qJ(4) + pkin(6);
t80 = cos(qJ(3));
t84 = pkin(3) * t80 + pkin(2);
t98 = t77 * t75 + t84 * t76 + pkin(1);
t97 = g(1) * t75;
t74 = qJ(3) + pkin(9);
t73 = qJ(5) + t74;
t68 = sin(t73);
t79 = sin(qJ(1));
t96 = t79 * t68;
t69 = cos(t73);
t95 = t79 * t69;
t71 = sin(t74);
t94 = t79 * t71;
t72 = cos(t74);
t93 = t79 * t72;
t78 = sin(qJ(3));
t92 = t79 * t78;
t91 = t79 * t80;
t81 = cos(qJ(1));
t90 = t81 * t68;
t89 = t81 * t69;
t88 = t81 * t71;
t87 = t81 * t72;
t86 = t81 * t78;
t85 = t81 * t80;
t82 = -g(2) * t79 + g(3) * t81;
t70 = -t78 * pkin(3) - qJ(2);
t67 = g(2) * t81 + g(3) * t79;
t1 = [0, t82, -t67, t82 * t76 - t97, t67, -g(1) * pkin(5) - g(2) * (t79 * pkin(1) - t81 * qJ(2)) - g(3) * (-t81 * pkin(1) - t79 * qJ(2)), 0, 0, 0, 0, 0, -t80 * t97 - g(2) * (t76 * t91 - t86) - g(3) * (-t76 * t85 - t92), t78 * t97 - g(2) * (-t76 * t92 - t85) - g(3) * (t76 * t86 - t91), -t72 * t97 - g(2) * (t76 * t93 - t88) - g(3) * (-t76 * t87 - t94), t71 * t97 - g(2) * (-t76 * t94 - t87) - g(3) * (t76 * t88 - t93), g(1) * t76 + t82 * t75, -g(1) * (t84 * t75 - t76 * t77 + pkin(5)) - g(2) * (t70 * t81 + t98 * t79) - g(3) * (t70 * t79 - t98 * t81), 0, 0, 0, 0, 0, -t69 * t97 - g(2) * (t76 * t95 - t90) - g(3) * (-t76 * t89 - t96), t68 * t97 - g(2) * (-t76 * t96 - t89) - g(3) * (t76 * t90 - t95);];
U_reg = t1;

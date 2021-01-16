% Calculate minimal parameter regressor of potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:34
% EndTime: 2021-01-15 12:56:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (86->40), mult. (112->59), div. (0->0), fcn. (118->8), ass. (0->28)
t68 = sin(pkin(8));
t85 = g(1) * t68;
t67 = qJ(3) + qJ(4);
t63 = sin(t67);
t71 = sin(qJ(1));
t84 = t71 * t63;
t64 = cos(t67);
t83 = t71 * t64;
t70 = sin(qJ(3));
t82 = t71 * t70;
t72 = cos(qJ(3));
t81 = t71 * t72;
t73 = cos(qJ(1));
t80 = t73 * t63;
t79 = t73 * t64;
t78 = t73 * t70;
t77 = t73 * t72;
t76 = t70 * pkin(3) + pkin(4) * t63 + qJ(2);
t75 = -g(2) * t71 + g(3) * t73;
t60 = t72 * pkin(3) + pkin(4) * t64 + pkin(2);
t66 = -qJ(5) - pkin(7) - pkin(6);
t69 = cos(pkin(8));
t74 = t60 * t69 - t66 * t68;
t65 = t71 * pkin(1);
t62 = g(2) * t73 + g(3) * t71;
t59 = -t64 * t85 - g(2) * (t69 * t83 - t80) - g(3) * (-t69 * t79 - t84);
t58 = t63 * t85 - g(2) * (-t69 * t84 - t79) - g(3) * (t69 * t80 - t83);
t1 = [0, t75, -t62, t75 * t69 - t85, t62, -g(1) * pkin(5) - g(2) * (-t73 * qJ(2) + t65) - g(3) * (-t73 * pkin(1) - t71 * qJ(2)), 0, 0, 0, 0, 0, -t72 * t85 - g(2) * (t69 * t81 - t78) - g(3) * (-t69 * t77 - t82), t70 * t85 - g(2) * (-t69 * t82 - t77) - g(3) * (t69 * t78 - t81), 0, 0, 0, 0, 0, t59, t58, t59, t58, g(1) * t69 + t75 * t68, -g(1) * (t68 * t60 + t69 * t66 + pkin(5)) - g(2) * t65 + (-g(2) * t74 + g(3) * t76) * t71 + (g(2) * t76 - g(3) * (-pkin(1) - t74)) * t73;];
U_reg = t1;

% Calculate minimal parameter regressor of potential energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:22
% EndTime: 2021-01-15 16:45:22
% DurationCPUTime: 0.23s
% Computational Cost: add. (120->67), mult. (281->109), div. (0->0), fcn. (348->10), ass. (0->43)
t74 = sin(pkin(5));
t76 = cos(pkin(5));
t79 = sin(qJ(3));
t80 = sin(qJ(2));
t82 = cos(qJ(3));
t86 = t80 * t82;
t68 = t74 * t86 + t76 * t79;
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t83 = cos(qJ(2));
t90 = t74 * t83;
t96 = g(3) * (-t68 * t78 - t81 * t90);
t95 = g(3) * (t68 * t81 - t78 * t90);
t73 = sin(pkin(9));
t75 = cos(pkin(9));
t87 = t76 * t83;
t64 = t73 * t80 - t75 * t87;
t59 = t64 * t78;
t66 = t73 * t87 + t75 * t80;
t61 = t66 * t78;
t94 = t73 * t76;
t93 = t74 * t79;
t92 = t74 * t80;
t91 = t74 * t82;
t89 = t75 * t76;
t88 = t76 * t80;
t85 = t82 * t83;
t63 = t73 * t88 - t75 * t83;
t51 = t63 * t79 + t73 * t91;
t65 = t73 * t83 + t75 * t88;
t52 = t65 * t79 + t75 * t91;
t67 = -t76 * t82 + t79 * t92;
t84 = g(1) * t51 - g(2) * t52 - g(3) * t67;
t77 = -qJ(5) - pkin(8);
t72 = pkin(4) * t81 + pkin(3);
t69 = t76 * t86 - t93;
t62 = t66 * t81;
t60 = t64 * t81;
t56 = -t73 * t69 + t75 * t85;
t55 = t75 * t69 + t73 * t85;
t54 = -t63 * t82 + t73 * t93;
t53 = t65 * t82 - t75 * t93;
t1 = [-g(3) * qJ(1), 0, g(1) * t63 - g(2) * t65 - g(3) * t92, g(1) * t66 + g(2) * t64 - g(3) * t90, 0, 0, 0, 0, 0, -g(1) * t54 - g(2) * t53 - g(3) * t68, -t84, 0, 0, 0, 0, 0, -g(1) * (t56 * t81 + t61) - g(2) * (t55 * t81 + t59) - t95, -g(1) * (-t56 * t78 + t62) - g(2) * (-t55 * t78 + t60) - t96, -g(1) * (t54 * t81 + t61) - g(2) * (t53 * t81 + t59) - t95, -g(1) * (-t54 * t78 + t62) - g(2) * (-t53 * t78 + t60) - t96, t84, -g(1) * (t54 * t72 + t51 * t77 + pkin(4) * t61 + (pkin(2) * t75 + pkin(7) * t94) * t83 + (-pkin(2) * t94 + pkin(7) * t75) * t80 + t75 * pkin(1)) - g(2) * (t53 * t72 - t52 * t77 + pkin(4) * t59 + (pkin(2) * t73 - pkin(7) * t89) * t83 + (pkin(2) * t89 + pkin(7) * t73) * t80 + t73 * pkin(1)) - g(3) * (t76 * pkin(6) - t67 * t77 + t68 * t72 + qJ(1)) + (-g(3) * (pkin(2) * t80 + (-pkin(4) * t78 - pkin(7)) * t83) + (-g(1) * t73 + g(2) * t75) * pkin(6)) * t74;];
U_reg = t1;

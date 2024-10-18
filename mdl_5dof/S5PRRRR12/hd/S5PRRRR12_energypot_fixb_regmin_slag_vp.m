% Calculate minimal parameter regressor of potential energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:19
% EndTime: 2024-09-28 18:08:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (157->63), mult. (145->109), div. (0->0), fcn. (164->22), ass. (0->45)
t68 = sin(pkin(5));
t87 = g(3) * t68;
t67 = sin(pkin(6));
t86 = t67 * t68;
t70 = cos(pkin(6));
t72 = sin(qJ(5));
t85 = t70 * t72;
t74 = cos(qJ(5));
t84 = t70 * t74;
t71 = cos(pkin(5));
t83 = t71 * t72;
t73 = sin(qJ(2));
t82 = t71 * t73;
t81 = t71 * t74;
t75 = cos(qJ(2));
t80 = t71 * t75;
t65 = qJ(2) + qJ(3);
t66 = sin(pkin(11));
t79 = t66 * t86;
t69 = cos(pkin(11));
t78 = t69 * t86;
t77 = t70 * t83;
t76 = t70 * t81;
t61 = pkin(5) - t65;
t60 = pkin(5) + t65;
t64 = qJ(4) + t65;
t63 = cos(t65);
t62 = sin(t65);
t59 = cos(t64);
t58 = sin(t64);
t57 = -qJ(4) + t61;
t56 = qJ(4) + t60;
t55 = cos(t60);
t54 = sin(t61);
t53 = cos(t56);
t52 = sin(t57);
t51 = cos(t61) / 0.2e1;
t50 = sin(t60) / 0.2e1;
t49 = cos(t57) / 0.2e1;
t48 = sin(t56) / 0.2e1;
t47 = t51 + t55 / 0.2e1;
t46 = t50 - t54 / 0.2e1;
t45 = t49 + t53 / 0.2e1;
t44 = t48 - t52 / 0.2e1;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t66 * t82 + t69 * t75) - g(2) * (t66 * t75 + t69 * t82) - t73 * t87, -g(1) * (-t66 * t80 - t69 * t73) - g(2) * (-t66 * t73 + t69 * t80) - t75 * t87, 0, -g(1) * (-t66 * t46 + t69 * t63) - g(2) * (t69 * t46 + t66 * t63) - g(3) * (t51 - t55 / 0.2e1), -g(1) * (-t66 * t47 - t69 * t62) - g(2) * (t69 * t47 - t66 * t62) - g(3) * (t50 + t54 / 0.2e1), 0, -g(1) * (-t66 * t44 + t69 * t59) - g(2) * (t69 * t44 + t66 * t59) - g(3) * (t49 - t53 / 0.2e1), -g(1) * (-t66 * t45 - t69 * t58) - g(2) * (t69 * t45 - t66 * t58) - g(3) * (t48 + t52 / 0.2e1), 0, 0, 0, 0, 0, -g(1) * ((-t66 * t77 + t69 * t74) * t59 + (-t66 * t81 - t69 * t85) * t58 + t72 * t79) - g(2) * ((t66 * t74 + t69 * t77) * t59 + (-t66 * t85 + t69 * t81) * t58 - t72 * t78) - g(3) * (t67 * t83 + (t58 * t74 + t59 * t85) * t68), -g(1) * ((-t66 * t76 - t69 * t72) * t59 + (t66 * t83 - t69 * t84) * t58 + t74 * t79) - g(2) * ((-t66 * t72 + t69 * t76) * t59 + (-t66 * t84 - t69 * t83) * t58 - t74 * t78) - g(3) * (t67 * t81 + (-t58 * t72 + t59 * t84) * t68);];
U_reg = t1;

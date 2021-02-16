% Calculate minimal parameter regressor of potential energy for
% S5RPRRP6
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:26
% EndTime: 2021-01-15 18:08:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (80->29), mult. (90->41), div. (0->0), fcn. (92->8), ass. (0->24)
t61 = sin(qJ(3));
t74 = g(3) * t61;
t60 = sin(qJ(4));
t64 = cos(qJ(3));
t73 = t60 * t64;
t63 = cos(qJ(4));
t72 = t63 * t64;
t71 = qJ(2) + pkin(5);
t70 = pkin(4) * t60 + pkin(6);
t58 = qJ(1) + pkin(8);
t56 = sin(t58);
t57 = cos(t58);
t69 = g(1) * t57 + g(2) * t56;
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t68 = -g(1) * t65 - g(2) * t62;
t55 = t63 * pkin(4) + pkin(3);
t59 = -qJ(5) - pkin(7);
t67 = t55 * t64 - t59 * t61 + pkin(2);
t66 = t68 * pkin(1);
t54 = -g(3) * t64 + t69 * t61;
t53 = -g(1) * (t56 * t60 + t57 * t72) - g(2) * (t56 * t72 - t57 * t60) - t63 * t74;
t52 = -g(1) * (t56 * t63 - t57 * t73) - g(2) * (-t56 * t73 - t57 * t63) + t60 * t74;
t1 = [0, t68, g(1) * t62 - g(2) * t65, -g(3) * t71 + t66, 0, 0, 0, 0, 0, -t69 * t64 - t74, t54, 0, 0, 0, 0, 0, t53, t52, t53, t52, -t54, -g(3) * (t61 * t55 + t64 * t59 + t71) + t66 + (-g(1) * t67 + g(2) * t70) * t57 + (-g(1) * t70 - g(2) * t67) * t56;];
U_reg = t1;

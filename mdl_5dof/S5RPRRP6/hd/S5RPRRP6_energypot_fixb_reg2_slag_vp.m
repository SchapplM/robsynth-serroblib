% Calculate inertial parameters regressor of potential energy for
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
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:15
% EndTime: 2019-12-31 18:43:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (128->49), mult. (130->62), div. (0->0), fcn. (125->8), ass. (0->29)
t64 = cos(qJ(4));
t52 = t64 * pkin(4) + pkin(3);
t59 = -qJ(5) - pkin(7);
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t80 = t52 * t65 - t59 * t62;
t60 = qJ(2) + pkin(5);
t79 = g(3) * t60;
t78 = g(3) * t62;
t58 = qJ(1) + pkin(8);
t53 = sin(t58);
t61 = sin(qJ(4));
t76 = t53 * t61;
t74 = t61 * t65;
t73 = t64 * t65;
t63 = sin(qJ(1));
t72 = t63 * pkin(1) + t53 * pkin(2);
t54 = cos(t58);
t66 = cos(qJ(1));
t71 = t66 * pkin(1) + t54 * pkin(2) + t53 * pkin(6);
t70 = -t54 * pkin(6) + t72;
t69 = pkin(3) * t65 + pkin(7) * t62;
t68 = g(1) * t54 + g(2) * t53;
t67 = -g(1) * t66 - g(2) * t63;
t48 = g(1) * t53 - g(2) * t54;
t47 = -g(3) * t65 + t68 * t62;
t46 = -g(1) * (t54 * t73 + t76) - g(2) * (t53 * t73 - t54 * t61) - t64 * t78;
t45 = -g(1) * (t53 * t64 - t54 * t74) - g(2) * (-t53 * t74 - t54 * t64) + t61 * t78;
t1 = [0, 0, 0, 0, 0, 0, t67, g(1) * t63 - g(2) * t66, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t68, t48, -g(3), t67 * pkin(1) - t79, 0, 0, 0, 0, 0, 0, -t68 * t65 - t78, t47, -t48, -g(1) * t71 - g(2) * t70 - t79, 0, 0, 0, 0, 0, 0, t46, t45, -t47, -g(1) * (t69 * t54 + t71) - g(2) * (t69 * t53 + t70) - g(3) * (t62 * pkin(3) - t65 * pkin(7) + t60), 0, 0, 0, 0, 0, 0, t46, t45, -t47, -g(1) * (pkin(4) * t76 + t71) - g(2) * (t80 * t53 + t72) - g(3) * (t62 * t52 + t65 * t59 + t60) + (-g(1) * t80 - g(2) * (-pkin(4) * t61 - pkin(6))) * t54;];
U_reg = t1;

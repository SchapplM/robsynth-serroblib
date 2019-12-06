% Calculate inertial parameters regressor of potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:24
% EndTime: 2019-12-05 17:36:24
% DurationCPUTime: 0.12s
% Computational Cost: add. (128->48), mult. (130->62), div. (0->0), fcn. (125->8), ass. (0->29)
t64 = cos(qJ(4));
t53 = t64 * pkin(4) + pkin(3);
t58 = sin(pkin(8));
t59 = cos(pkin(8));
t60 = -qJ(5) - pkin(6);
t79 = t53 * t59 - t58 * t60;
t78 = g(1) * t58;
t61 = qJ(2) + pkin(5);
t77 = g(1) * t61;
t57 = qJ(1) + pkin(7);
t54 = sin(t57);
t76 = g(2) * t54;
t55 = cos(t57);
t62 = sin(qJ(4));
t74 = t55 * t62;
t72 = t59 * t62;
t71 = t59 * t64;
t65 = cos(qJ(1));
t70 = t65 * pkin(1) + t55 * pkin(2) + t54 * qJ(3);
t63 = sin(qJ(1));
t69 = -t63 * pkin(1) + t55 * qJ(3);
t68 = pkin(3) * t59 + pkin(6) * t58;
t67 = -g(3) * t55 + t76;
t66 = g(2) * t63 - g(3) * t65;
t49 = g(2) * t55 + g(3) * t54;
t48 = g(1) * t59 + t67 * t58;
t47 = -t64 * t78 - g(2) * (-t54 * t71 + t74) - g(3) * (t54 * t62 + t55 * t71);
t46 = t62 * t78 - g(2) * (t54 * t72 + t55 * t64) - g(3) * (t54 * t64 - t55 * t72);
t1 = [0, 0, 0, 0, 0, 0, t66, g(2) * t65 + g(3) * t63, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t67, t49, -g(1), t66 * pkin(1) - t77, 0, 0, 0, 0, 0, 0, t67 * t59 - t78, -t48, -t49, -t77 - g(2) * (-t54 * pkin(2) + t69) - g(3) * t70, 0, 0, 0, 0, 0, 0, t47, t46, t48, -g(1) * (t58 * pkin(3) - t59 * pkin(6) + t61) - g(2) * t69 - g(3) * (t68 * t55 + t70) - (-pkin(2) - t68) * t76, 0, 0, 0, 0, 0, 0, t47, t46, t48, -g(1) * (t58 * t53 + t59 * t60 + t61) - g(2) * (pkin(4) * t74 + t69) - g(3) * (t79 * t55 + t70) + (-g(2) * (-pkin(2) - t79) - g(3) * pkin(4) * t62) * t54;];
U_reg = t1;

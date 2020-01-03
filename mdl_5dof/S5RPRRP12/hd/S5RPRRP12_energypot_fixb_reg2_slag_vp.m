% Calculate inertial parameters regressor of potential energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:19
% EndTime: 2019-12-31 18:57:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (77->49), mult. (132->57), div. (0->0), fcn. (127->6), ass. (0->30)
t54 = cos(qJ(4));
t44 = t54 * pkin(4) + pkin(3);
t50 = -qJ(5) - pkin(7);
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t72 = t44 * t52 + t50 * t55;
t71 = g(3) * pkin(5);
t70 = pkin(2) + pkin(5);
t56 = cos(qJ(1));
t69 = g(2) * t56;
t68 = g(3) * t55;
t51 = sin(qJ(4));
t53 = sin(qJ(1));
t65 = t53 * t51;
t64 = t53 * t54;
t63 = t56 * t51;
t62 = t56 * t54;
t46 = t53 * pkin(6);
t47 = t53 * pkin(1);
t61 = t46 + t47;
t60 = t56 * pkin(1) + t53 * qJ(2);
t59 = t56 * pkin(6) + t60;
t58 = -t56 * qJ(2) + t47;
t57 = pkin(3) * t52 - pkin(7) * t55;
t42 = g(1) * t53 - t69;
t43 = g(1) * t56 + g(2) * t53;
t41 = -g(3) * t52 + t42 * t55;
t40 = -g(1) * (t52 * t64 + t63) - g(2) * (-t52 * t62 + t65) - t54 * t68;
t39 = -g(1) * (-t52 * t65 + t62) - g(2) * (t52 * t63 + t64) + t51 * t68;
t1 = [0, 0, 0, 0, 0, 0, -t43, t42, -g(3), -t71, 0, 0, 0, 0, 0, 0, -g(3), t43, -t42, -g(1) * t60 - g(2) * t58 - t71, 0, 0, 0, 0, 0, 0, -t42 * t52 - t68, -t41, -t43, -g(1) * t59 - g(2) * (t46 + t58) - g(3) * t70, 0, 0, 0, 0, 0, 0, t40, t39, t41, -g(1) * (t57 * t53 + t59) - g(2) * t61 - g(3) * (t55 * pkin(3) + t52 * pkin(7) + t70) - (-qJ(2) - t57) * t69, 0, 0, 0, 0, 0, 0, t40, t39, t41, -g(1) * (t72 * t53 + t59) - g(2) * (pkin(4) * t65 + t61) - g(3) * (t55 * t44 - t52 * t50 + t70) + (-g(1) * pkin(4) * t51 - g(2) * (-qJ(2) - t72)) * t56;];
U_reg = t1;

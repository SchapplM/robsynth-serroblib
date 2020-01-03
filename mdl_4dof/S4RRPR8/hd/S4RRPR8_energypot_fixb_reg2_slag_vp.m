% Calculate inertial parameters regressor of potential energy for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:22
% EndTime: 2019-12-31 17:08:22
% DurationCPUTime: 0.07s
% Computational Cost: add. (54->35), mult. (105->43), div. (0->0), fcn. (101->6), ass. (0->23)
t53 = sin(qJ(1));
t52 = sin(qJ(2));
t63 = qJ(3) * t52;
t55 = cos(qJ(2));
t66 = t53 * t55;
t70 = pkin(2) * t66 + t53 * t63;
t56 = cos(qJ(1));
t59 = g(1) * t56 + g(2) * t53;
t69 = g(3) * pkin(4);
t65 = t55 * t56;
t64 = t56 * pkin(1) + t53 * pkin(5);
t48 = t53 * pkin(1);
t62 = -t56 * pkin(5) + t48;
t61 = pkin(2) * t65 + t56 * t63 + t64;
t60 = t52 * pkin(2) - t55 * qJ(3) + pkin(4);
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t58 = t55 * t51 - t52 * t54;
t57 = t52 * t51 + t55 * t54;
t41 = g(1) * t53 - g(2) * t56;
t40 = -g(3) * t52 - t59 * t55;
t39 = -g(3) * t55 + t59 * t52;
t1 = [0, 0, 0, 0, 0, 0, -t59, t41, -g(3), -t69, 0, 0, 0, 0, 0, 0, t40, t39, -t41, -g(1) * t64 - g(2) * t62 - t69, 0, 0, 0, 0, 0, 0, t40, -t41, -t39, -g(1) * t61 - g(2) * (t62 + t70) - g(3) * t60, 0, 0, 0, 0, 0, 0, g(3) * t58 - t59 * t57, g(3) * t57 + t59 * t58, t41, -g(1) * (pkin(3) * t65 - t53 * pkin(6) + t61) - g(2) * (pkin(3) * t66 + t48 + (-pkin(5) + pkin(6)) * t56 + t70) - g(3) * (t52 * pkin(3) + t60);];
U_reg = t1;

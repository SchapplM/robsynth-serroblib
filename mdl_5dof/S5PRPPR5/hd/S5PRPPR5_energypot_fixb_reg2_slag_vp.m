% Calculate inertial parameters regressor of potential energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:36
% EndTime: 2019-12-31 17:38:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (100->52), mult. (204->72), div. (0->0), fcn. (216->8), ass. (0->33)
t49 = sin(pkin(8));
t51 = cos(pkin(8));
t54 = sin(qJ(2));
t56 = cos(qJ(2));
t74 = t56 * t49 - t54 * t51;
t73 = g(3) * t74;
t72 = g(3) * qJ(1);
t50 = sin(pkin(7));
t71 = t50 * t56;
t52 = cos(pkin(7));
t70 = t52 * t56;
t67 = t52 * pkin(1) + t50 * pkin(5);
t66 = qJ(3) * t54;
t65 = t50 * pkin(1) - t52 * pkin(5);
t64 = pkin(2) * t70 + t52 * t66 + t67;
t63 = g(1) * t52 + g(2) * t50;
t62 = t54 * pkin(2) - t56 * qJ(3) + qJ(1);
t31 = t54 * t49 + t56 * t51;
t61 = pkin(2) * t71 + t50 * t66 + t65;
t60 = t54 * pkin(3) + t62;
t27 = t74 * t50;
t29 = t74 * t52;
t59 = g(1) * t29 + g(2) * t27 + g(3) * t31;
t58 = pkin(3) * t71 + t52 * qJ(4) + t61;
t57 = pkin(3) * t70 - t50 * qJ(4) + t64;
t55 = cos(qJ(5));
t53 = sin(qJ(5));
t33 = g(1) * t50 - g(2) * t52;
t30 = t31 * t52;
t28 = t31 * t50;
t26 = -g(3) * t54 - t63 * t56;
t25 = -g(3) * t56 + t63 * t54;
t1 = [0, 0, 0, 0, 0, 0, -t63, t33, -g(3), -t72, 0, 0, 0, 0, 0, 0, t26, t25, -t33, -g(1) * t67 - g(2) * t65 - t72, 0, 0, 0, 0, 0, 0, t26, -t33, -t25, -g(1) * t64 - g(2) * t61 - g(3) * t62, 0, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t28 + t73, t59, t33, -g(1) * t57 - g(2) * t58 - g(3) * t60, 0, 0, 0, 0, 0, 0, -g(1) * (t30 * t55 - t50 * t53) - g(2) * (t28 * t55 + t52 * t53) + t55 * t73, -g(1) * (-t30 * t53 - t50 * t55) - g(2) * (-t28 * t53 + t52 * t55) - t53 * t73, -t59, -g(1) * (t30 * pkin(4) + t29 * pkin(6) + t57) - g(2) * (t28 * pkin(4) + t27 * pkin(6) + t58) - g(3) * (-pkin(4) * t74 + t31 * pkin(6) + t60);];
U_reg = t1;

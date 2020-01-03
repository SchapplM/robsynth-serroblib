% Calculate inertial parameters regressor of potential energy for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:39
% EndTime: 2019-12-31 20:15:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (99->40), mult. (78->41), div. (0->0), fcn. (65->8), ass. (0->22)
t57 = pkin(6) + pkin(5);
t52 = sin(qJ(4));
t64 = pkin(4) * t52;
t63 = g(3) * t57;
t51 = qJ(1) + qJ(2);
t45 = sin(t51);
t53 = sin(qJ(1));
t62 = t53 * pkin(1) + t45 * pkin(2);
t61 = pkin(3) + t57;
t47 = cos(t51);
t55 = cos(qJ(1));
t60 = t55 * pkin(1) + t47 * pkin(2) + t45 * qJ(3);
t59 = -t47 * qJ(3) + t62;
t39 = g(1) * t45 - g(2) * t47;
t58 = -g(1) * t55 - g(2) * t53;
t56 = -pkin(8) - pkin(7);
t54 = cos(qJ(4));
t50 = qJ(4) + qJ(5);
t46 = cos(t50);
t44 = sin(t50);
t40 = g(1) * t47 + g(2) * t45;
t1 = [0, 0, 0, 0, 0, 0, t58, g(1) * t53 - g(2) * t55, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t40, t39, -g(3), t58 * pkin(1) - t63, 0, 0, 0, 0, 0, 0, -g(3), t40, -t39, -g(1) * t60 - g(2) * t59 - t63, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t39 * t52, g(3) * t52 - t39 * t54, -t40, -g(1) * (t47 * pkin(7) + t60) - g(2) * (t45 * pkin(7) + t59) - g(3) * t61, 0, 0, 0, 0, 0, 0, -g(3) * t46 - t39 * t44, g(3) * t44 - t39 * t46, -t40, -g(1) * (t45 * t64 - t47 * t56 + t60) - g(2) * (-t45 * t56 + (-qJ(3) - t64) * t47 + t62) - g(3) * (t54 * pkin(4) + t61);];
U_reg = t1;

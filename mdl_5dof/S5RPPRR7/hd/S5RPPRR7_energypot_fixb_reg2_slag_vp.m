% Calculate inertial parameters regressor of potential energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:52
% EndTime: 2019-12-31 17:59:52
% DurationCPUTime: 0.07s
% Computational Cost: add. (105->45), mult. (96->53), div. (0->0), fcn. (87->8), ass. (0->27)
t49 = qJ(1) + pkin(8);
t46 = cos(t49);
t68 = g(2) * t46;
t50 = qJ(2) + pkin(5);
t67 = g(3) * t50;
t55 = cos(qJ(4));
t66 = g(3) * t55;
t51 = sin(qJ(5));
t52 = sin(qJ(4));
t65 = t51 * t52;
t54 = cos(qJ(5));
t64 = t52 * t54;
t45 = sin(t49);
t53 = sin(qJ(1));
t63 = t53 * pkin(1) + t45 * pkin(2);
t62 = pkin(3) + t50;
t56 = cos(qJ(1));
t61 = t56 * pkin(1) + t46 * pkin(2) + t45 * qJ(3);
t60 = t46 * pkin(6) + t61;
t59 = -t46 * qJ(3) + t63;
t58 = pkin(4) * t52 - pkin(7) * t55;
t38 = g(1) * t45 - t68;
t57 = -g(1) * t56 - g(2) * t53;
t41 = t45 * pkin(6);
t39 = g(1) * t46 + g(2) * t45;
t37 = -g(3) * t52 + t38 * t55;
t1 = [0, 0, 0, 0, 0, 0, t57, g(1) * t53 - g(2) * t56, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t39, t38, -g(3), t57 * pkin(1) - t67, 0, 0, 0, 0, 0, 0, -g(3), t39, -t38, -g(1) * t61 - g(2) * t59 - t67, 0, 0, 0, 0, 0, 0, -t38 * t52 - t66, -t37, -t39, -g(1) * t60 - g(2) * (t41 + t59) - g(3) * t62, 0, 0, 0, 0, 0, 0, -g(1) * (t45 * t64 + t46 * t51) - g(2) * (t45 * t51 - t46 * t64) - t54 * t66, -g(1) * (-t45 * t65 + t46 * t54) - g(2) * (t45 * t54 + t46 * t65) + t51 * t66, t37, -g(1) * (t58 * t45 + t60) - g(2) * (t41 + t63) - g(3) * (t55 * pkin(4) + t52 * pkin(7) + t62) - (-qJ(3) - t58) * t68;];
U_reg = t1;

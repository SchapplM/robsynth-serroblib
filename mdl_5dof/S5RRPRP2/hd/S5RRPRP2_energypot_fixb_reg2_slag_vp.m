% Calculate inertial parameters regressor of potential energy for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:42
% EndTime: 2019-12-31 19:49:42
% DurationCPUTime: 0.06s
% Computational Cost: add. (121->38), mult. (81->41), div. (0->0), fcn. (68->8), ass. (0->24)
t69 = pkin(6) + pkin(5);
t55 = qJ(3) + t69;
t68 = g(3) * t55;
t56 = qJ(1) + qJ(2);
t51 = sin(t56);
t58 = sin(qJ(1));
t67 = t58 * pkin(1) + pkin(2) * t51;
t52 = cos(t56);
t60 = cos(qJ(1));
t66 = t60 * pkin(1) + pkin(2) * t52;
t50 = pkin(8) + t56;
t46 = sin(t50);
t47 = cos(t50);
t65 = t47 * pkin(3) + t46 * pkin(7) + t66;
t64 = g(1) * t47 + g(2) * t46;
t63 = -g(1) * t60 - g(2) * t58;
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t62 = pkin(4) * t59 + qJ(5) * t57;
t61 = t46 * pkin(3) - t47 * pkin(7) + t67;
t41 = g(1) * t46 - g(2) * t47;
t40 = -g(3) * t57 - t64 * t59;
t39 = -g(3) * t59 + t64 * t57;
t1 = [0, 0, 0, 0, 0, 0, t63, g(1) * t58 - g(2) * t60, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t52 - g(2) * t51, g(1) * t51 - g(2) * t52, -g(3), t63 * pkin(1) - g(3) * t69, 0, 0, 0, 0, 0, 0, -t64, t41, -g(3), -g(1) * t66 - g(2) * t67 - t68, 0, 0, 0, 0, 0, 0, t40, t39, -t41, -g(1) * t65 - g(2) * t61 - t68, 0, 0, 0, 0, 0, 0, t40, -t41, -t39, -g(1) * (t62 * t47 + t65) - g(2) * (t62 * t46 + t61) - g(3) * (t57 * pkin(4) - t59 * qJ(5) + t55);];
U_reg = t1;

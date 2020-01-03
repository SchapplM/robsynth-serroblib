% Calculate inertial parameters regressor of potential energy for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:56
% EndTime: 2019-12-31 18:40:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (121->38), mult. (81->41), div. (0->0), fcn. (68->8), ass. (0->24)
t69 = qJ(2) + pkin(5);
t56 = pkin(6) + t69;
t70 = g(3) * t56;
t57 = qJ(1) + pkin(8);
t51 = sin(t57);
t59 = sin(qJ(1));
t68 = t59 * pkin(1) + pkin(2) * t51;
t52 = cos(t57);
t61 = cos(qJ(1));
t67 = t61 * pkin(1) + pkin(2) * t52;
t53 = qJ(3) + t57;
t49 = sin(t53);
t50 = cos(t53);
t66 = t50 * pkin(3) + t49 * pkin(7) + t67;
t65 = g(1) * t50 + g(2) * t49;
t64 = -g(1) * t61 - g(2) * t59;
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t63 = pkin(4) * t60 + qJ(5) * t58;
t62 = t49 * pkin(3) - t50 * pkin(7) + t68;
t42 = g(1) * t49 - g(2) * t50;
t41 = -g(3) * t58 - t65 * t60;
t40 = -g(3) * t60 + t65 * t58;
t1 = [0, 0, 0, 0, 0, 0, t64, g(1) * t59 - g(2) * t61, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t52 - g(2) * t51, g(1) * t51 - g(2) * t52, -g(3), t64 * pkin(1) - g(3) * t69, 0, 0, 0, 0, 0, 0, -t65, t42, -g(3), -g(1) * t67 - g(2) * t68 - t70, 0, 0, 0, 0, 0, 0, t41, t40, -t42, -g(1) * t66 - g(2) * t62 - t70, 0, 0, 0, 0, 0, 0, t41, -t42, -t40, -g(1) * (t63 * t50 + t66) - g(2) * (t63 * t49 + t62) - g(3) * (t58 * pkin(4) - t60 * qJ(5) + t56);];
U_reg = t1;

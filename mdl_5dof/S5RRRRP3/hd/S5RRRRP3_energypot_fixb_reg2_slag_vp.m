% Calculate inertial parameters regressor of potential energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:24
% EndTime: 2019-12-31 21:49:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (121->38), mult. (81->41), div. (0->0), fcn. (68->8), ass. (0->24)
t71 = pkin(6) + pkin(5);
t57 = pkin(7) + t71;
t70 = g(3) * t57;
t58 = qJ(1) + qJ(2);
t52 = sin(t58);
t60 = sin(qJ(1));
t69 = t60 * pkin(1) + pkin(2) * t52;
t53 = cos(t58);
t62 = cos(qJ(1));
t68 = t62 * pkin(1) + pkin(2) * t53;
t54 = qJ(3) + t58;
t50 = sin(t54);
t51 = cos(t54);
t67 = t51 * pkin(3) + t50 * pkin(8) + t68;
t66 = g(1) * t51 + g(2) * t50;
t65 = -g(1) * t62 - g(2) * t60;
t59 = sin(qJ(4));
t61 = cos(qJ(4));
t64 = pkin(4) * t61 + qJ(5) * t59;
t63 = t50 * pkin(3) - t51 * pkin(8) + t69;
t43 = g(1) * t50 - g(2) * t51;
t42 = -g(3) * t59 - t66 * t61;
t41 = -g(3) * t61 + t66 * t59;
t1 = [0, 0, 0, 0, 0, 0, t65, g(1) * t60 - g(2) * t62, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t52, g(1) * t52 - g(2) * t53, -g(3), t65 * pkin(1) - g(3) * t71, 0, 0, 0, 0, 0, 0, -t66, t43, -g(3), -g(1) * t68 - g(2) * t69 - t70, 0, 0, 0, 0, 0, 0, t42, t41, -t43, -g(1) * t67 - g(2) * t63 - t70, 0, 0, 0, 0, 0, 0, t42, -t43, -t41, -g(1) * (t64 * t51 + t67) - g(2) * (t64 * t50 + t63) - g(3) * (t59 * pkin(4) - t61 * qJ(5) + t57);];
U_reg = t1;

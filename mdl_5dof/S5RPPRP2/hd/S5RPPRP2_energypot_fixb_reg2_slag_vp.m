% Calculate inertial parameters regressor of potential energy for
% S5RPPRP2
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:23
% EndTime: 2019-12-31 17:49:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->42), mult. (93->45), div. (0->0), fcn. (80->8), ass. (0->26)
t60 = qJ(2) + pkin(5);
t70 = g(3) * t60;
t59 = cos(pkin(8));
t48 = t59 * pkin(3) + pkin(2);
t57 = qJ(1) + pkin(7);
t50 = sin(t57);
t52 = cos(t57);
t62 = sin(qJ(1));
t54 = t62 * pkin(1);
t61 = -pkin(6) - qJ(3);
t69 = t50 * t48 + t52 * t61 + t54;
t58 = sin(pkin(8));
t68 = t58 * pkin(3) + t60;
t63 = cos(qJ(1));
t55 = t63 * pkin(1);
t67 = t52 * t48 - t50 * t61 + t55;
t66 = g(1) * t52 + g(2) * t50;
t65 = -g(1) * t63 - g(2) * t62;
t56 = pkin(8) + qJ(4);
t49 = sin(t56);
t51 = cos(t56);
t64 = pkin(4) * t51 + qJ(5) * t49;
t43 = g(1) * t50 - g(2) * t52;
t42 = -g(3) * t49 - t66 * t51;
t41 = -g(3) * t51 + t66 * t49;
t1 = [0, 0, 0, 0, 0, 0, t65, g(1) * t62 - g(2) * t63, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t66, t43, -g(3), t65 * pkin(1) - t70, 0, 0, 0, 0, 0, 0, -g(3) * t58 - t66 * t59, -g(3) * t59 + t66 * t58, -t43, -g(1) * (t52 * pkin(2) + t50 * qJ(3) + t55) - g(2) * (t50 * pkin(2) - t52 * qJ(3) + t54) - t70, 0, 0, 0, 0, 0, 0, t42, t41, -t43, -g(1) * t67 - g(2) * t69 - g(3) * t68, 0, 0, 0, 0, 0, 0, t42, -t43, -t41, -g(1) * (t64 * t52 + t67) - g(2) * (t64 * t50 + t69) - g(3) * (t49 * pkin(4) - t51 * qJ(5) + t68);];
U_reg = t1;

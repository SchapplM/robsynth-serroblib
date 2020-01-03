% Calculate inertial parameters regressor of potential energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:05
% EndTime: 2019-12-31 19:06:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (87->41), mult. (130->47), div. (0->0), fcn. (141->8), ass. (0->25)
t73 = g(3) * pkin(5);
t62 = -pkin(6) + pkin(5);
t72 = g(3) * t62;
t71 = cos(qJ(3));
t70 = sin(qJ(1));
t69 = sin(qJ(3));
t60 = cos(qJ(1));
t68 = t60 * pkin(1) + t70 * qJ(2);
t67 = t60 * pkin(2) + t68;
t66 = t70 * pkin(1) - t60 * qJ(2);
t65 = t70 * pkin(2) + t66;
t42 = -t60 * t71 - t70 * t69;
t43 = t60 * t69 - t70 * t71;
t64 = g(1) * t43 - g(2) * t42;
t63 = g(1) * t42 + g(2) * t43;
t61 = -pkin(8) - pkin(7);
t59 = cos(qJ(4));
t58 = sin(qJ(4));
t57 = qJ(4) + qJ(5);
t51 = cos(t57);
t50 = sin(t57);
t49 = t59 * pkin(4) + pkin(3);
t45 = -g(1) * t60 - g(2) * t70;
t44 = g(1) * t70 - g(2) * t60;
t1 = [0, 0, 0, 0, 0, 0, t45, t44, -g(3), -t73, 0, 0, 0, 0, 0, 0, t45, -g(3), -t44, -g(1) * t68 - g(2) * t66 - t73, 0, 0, 0, 0, 0, 0, t63, t64, g(3), -g(1) * t67 - g(2) * t65 - t72, 0, 0, 0, 0, 0, 0, g(3) * t58 + t63 * t59, g(3) * t59 - t63 * t58, -t64, -g(1) * (-t42 * pkin(3) + t43 * pkin(7) + t67) - g(2) * (-t43 * pkin(3) - t42 * pkin(7) + t65) - t72, 0, 0, 0, 0, 0, 0, g(3) * t50 + t63 * t51, g(3) * t51 - t63 * t50, -t64, -g(1) * (-t42 * t49 - t43 * t61 + t67) - g(2) * (t42 * t61 - t43 * t49 + t65) - g(3) * (-t58 * pkin(4) + t62);];
U_reg = t1;

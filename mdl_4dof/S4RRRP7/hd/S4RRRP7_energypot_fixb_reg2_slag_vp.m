% Calculate inertial parameters regressor of potential energy for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:03
% EndTime: 2019-12-31 17:21:03
% DurationCPUTime: 0.07s
% Computational Cost: add. (63->39), mult. (128->50), div. (0->0), fcn. (130->6), ass. (0->27)
t74 = g(3) * pkin(4);
t57 = sin(qJ(2));
t73 = pkin(6) * t57;
t72 = g(3) * t57;
t58 = sin(qJ(1));
t60 = cos(qJ(2));
t71 = t58 * t60;
t56 = sin(qJ(3));
t61 = cos(qJ(1));
t70 = t61 * t56;
t59 = cos(qJ(3));
t69 = t61 * t59;
t68 = t61 * pkin(1) + t58 * pkin(5);
t67 = t58 * pkin(1) - t61 * pkin(5);
t66 = t68 + (pkin(2) * t60 + t73) * t61;
t65 = t57 * pkin(2) - t60 * pkin(6) + pkin(4);
t64 = g(1) * t61 + g(2) * t58;
t63 = pkin(2) * t71 + t58 * t73 + t67;
t41 = t56 * t71 + t69;
t43 = -t58 * t59 + t60 * t70;
t62 = g(1) * t43 + g(2) * t41 + t56 * t72;
t45 = g(1) * t58 - g(2) * t61;
t44 = t58 * t56 + t60 * t69;
t42 = t59 * t71 - t70;
t40 = -g(3) * t60 + t64 * t57;
t39 = -g(1) * t44 - g(2) * t42 - t59 * t72;
t1 = [0, 0, 0, 0, 0, 0, -t64, t45, -g(3), -t74, 0, 0, 0, 0, 0, 0, -t64 * t60 - t72, t40, -t45, -g(1) * t68 - g(2) * t67 - t74, 0, 0, 0, 0, 0, 0, t39, t62, -t40, -g(1) * t66 - g(2) * t63 - g(3) * t65, 0, 0, 0, 0, 0, 0, t39, -t40, -t62, -g(1) * (t44 * pkin(3) + t43 * qJ(4) + t66) - g(2) * (t42 * pkin(3) + t41 * qJ(4) + t63) - g(3) * ((pkin(3) * t59 + qJ(4) * t56) * t57 + t65);];
U_reg = t1;

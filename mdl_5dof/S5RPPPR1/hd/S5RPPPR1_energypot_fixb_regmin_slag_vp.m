% Calculate minimal parameter regressor of potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:43
% EndTime: 2022-01-20 09:12:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (86->38), mult. (84->61), div. (0->0), fcn. (83->10), ass. (0->24)
t63 = sin(pkin(8));
t78 = g(3) * t63;
t66 = qJ(2) + pkin(5);
t77 = g(3) * t66;
t61 = qJ(1) + pkin(7);
t55 = sin(t61);
t65 = cos(pkin(8));
t76 = t55 * t65;
t57 = cos(t61);
t75 = t57 * t65;
t62 = sin(pkin(9));
t74 = t62 * t65;
t64 = cos(pkin(9));
t73 = t64 * t65;
t68 = cos(qJ(1));
t72 = t68 * pkin(1) + t57 * pkin(2) + t55 * qJ(3);
t67 = sin(qJ(1));
t71 = t67 * pkin(1) + t55 * pkin(2) - t57 * qJ(3);
t70 = -g(1) * t68 - g(2) * t67;
t69 = pkin(3) * t65 + qJ(4) * t63;
t60 = pkin(9) + qJ(5);
t56 = cos(t60);
t54 = sin(t60);
t1 = [0, t70, g(1) * t67 - g(2) * t68, t70 * pkin(1) - t77, -t78 + (-g(1) * t57 - g(2) * t55) * t65, -g(1) * t55 + g(2) * t57, -g(1) * t72 - g(2) * t71 - t77, -g(1) * (t55 * t62 + t57 * t73) - g(2) * (t55 * t73 - t57 * t62) - t64 * t78, -g(1) * (t55 * t64 - t57 * t74) - g(2) * (-t55 * t74 - t57 * t64) + t62 * t78, -g(1) * (t69 * t57 + t72) - g(2) * (t69 * t55 + t71) - g(3) * (t63 * pkin(3) - t65 * qJ(4) + t66), 0, 0, 0, 0, 0, -g(1) * (t55 * t54 + t56 * t75) - g(2) * (-t57 * t54 + t56 * t76) - t56 * t78, -g(1) * (-t54 * t75 + t55 * t56) - g(2) * (-t54 * t76 - t57 * t56) + t54 * t78;];
U_reg = t1;

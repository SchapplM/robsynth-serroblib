% Calculate inertial parameters regressor of potential energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:53
% EndTime: 2019-12-31 18:12:53
% DurationCPUTime: 0.09s
% Computational Cost: add. (105->45), mult. (117->46), div. (0->0), fcn. (104->6), ass. (0->24)
t60 = pkin(7) + qJ(3);
t56 = sin(t60);
t57 = cos(t60);
t77 = pkin(3) * t57 + qJ(4) * t56;
t65 = cos(qJ(1));
t76 = t77 * t65;
t75 = g(3) * pkin(5);
t61 = sin(pkin(7));
t73 = t61 * pkin(2) + pkin(5);
t62 = cos(pkin(7));
t54 = t62 * pkin(2) + pkin(1);
t63 = -pkin(6) - qJ(2);
t64 = sin(qJ(1));
t72 = t64 * t54 + t65 * t63;
t70 = qJ(5) * t57;
t48 = t65 * t54;
t69 = -t64 * t63 + t48;
t68 = t77 * t64 + t72;
t67 = g(1) * t65 + g(2) * t64;
t66 = t56 * pkin(3) - t57 * qJ(4) + t73;
t49 = g(1) * t64 - g(2) * t65;
t44 = g(3) * t56 + t67 * t57;
t43 = -g(3) * t57 + t67 * t56;
t1 = [0, 0, 0, 0, 0, 0, -t67, t49, -g(3), -t75, 0, 0, 0, 0, 0, 0, -g(3) * t61 - t67 * t62, -g(3) * t62 + t67 * t61, -t49, -g(1) * (t65 * pkin(1) + t64 * qJ(2)) - g(2) * (t64 * pkin(1) - t65 * qJ(2)) - t75, 0, 0, 0, 0, 0, 0, -t44, t43, -t49, -g(1) * t69 - g(2) * t72 - g(3) * t73, 0, 0, 0, 0, 0, 0, -t49, t44, -t43, -g(1) * (t69 + t76) - g(2) * t68 - g(3) * t66, 0, 0, 0, 0, 0, 0, -t49, -t43, -t44, -g(1) * (t65 * t70 + t48 + (pkin(4) - t63) * t64 + t76) - g(2) * (-t65 * pkin(4) + t64 * t70 + t68) - g(3) * (t56 * qJ(5) + t66);];
U_reg = t1;

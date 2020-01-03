% Calculate inertial parameters regressor of potential energy for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:12
% EndTime: 2019-12-31 21:51:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->42), mult. (93->45), div. (0->0), fcn. (80->8), ass. (0->26)
t69 = pkin(6) + pkin(5);
t76 = g(3) * t69;
t64 = sin(qJ(3));
t75 = t64 * pkin(3) + t69;
t66 = cos(qJ(3));
t54 = t66 * pkin(3) + pkin(2);
t63 = qJ(1) + qJ(2);
t56 = sin(t63);
t58 = cos(t63);
t65 = sin(qJ(1));
t60 = t65 * pkin(1);
t68 = -pkin(8) - pkin(7);
t74 = t56 * t54 + t58 * t68 + t60;
t67 = cos(qJ(1));
t61 = t67 * pkin(1);
t73 = t58 * t54 - t56 * t68 + t61;
t72 = g(1) * t58 + g(2) * t56;
t71 = -g(1) * t67 - g(2) * t65;
t62 = qJ(3) + qJ(4);
t55 = sin(t62);
t57 = cos(t62);
t70 = pkin(4) * t57 + qJ(5) * t55;
t49 = g(1) * t56 - g(2) * t58;
t48 = -g(3) * t55 - t72 * t57;
t47 = -g(3) * t57 + t72 * t55;
t1 = [0, 0, 0, 0, 0, 0, t71, g(1) * t65 - g(2) * t67, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t72, t49, -g(3), t71 * pkin(1) - t76, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t72 * t66, -g(3) * t66 + t72 * t64, -t49, -g(1) * (t58 * pkin(2) + t56 * pkin(7) + t61) - g(2) * (t56 * pkin(2) - t58 * pkin(7) + t60) - t76, 0, 0, 0, 0, 0, 0, t48, t47, -t49, -g(1) * t73 - g(2) * t74 - g(3) * t75, 0, 0, 0, 0, 0, 0, t48, -t49, -t47, -g(1) * (t70 * t58 + t73) - g(2) * (t70 * t56 + t74) - g(3) * (t55 * pkin(4) - t57 * qJ(5) + t75);];
U_reg = t1;

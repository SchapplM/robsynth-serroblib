% Calculate inertial parameters regressor of potential energy for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:49:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->42), mult. (93->45), div. (0->0), fcn. (80->8), ass. (0->26)
t68 = pkin(6) + pkin(5);
t75 = g(3) * t68;
t64 = sin(qJ(3));
t74 = t64 * pkin(3) + t68;
t66 = cos(qJ(3));
t53 = t66 * pkin(3) + pkin(2);
t62 = qJ(1) + qJ(2);
t56 = sin(t62);
t57 = cos(t62);
t65 = sin(qJ(1));
t59 = t65 * pkin(1);
t63 = -qJ(4) - pkin(7);
t73 = t56 * t53 + t57 * t63 + t59;
t67 = cos(qJ(1));
t60 = t67 * pkin(1);
t72 = t57 * t53 - t56 * t63 + t60;
t71 = g(1) * t57 + g(2) * t56;
t70 = -g(1) * t67 - g(2) * t65;
t61 = qJ(3) + pkin(8);
t54 = sin(t61);
t55 = cos(t61);
t69 = pkin(4) * t55 + qJ(5) * t54;
t48 = g(1) * t56 - g(2) * t57;
t47 = -g(3) * t54 - t71 * t55;
t46 = -g(3) * t55 + t71 * t54;
t1 = [0, 0, 0, 0, 0, 0, t70, g(1) * t65 - g(2) * t67, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t71, t48, -g(3), t70 * pkin(1) - t75, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t71 * t66, -g(3) * t66 + t71 * t64, -t48, -g(1) * (t57 * pkin(2) + t56 * pkin(7) + t60) - g(2) * (t56 * pkin(2) - t57 * pkin(7) + t59) - t75, 0, 0, 0, 0, 0, 0, t47, t46, -t48, -g(1) * t72 - g(2) * t73 - g(3) * t74, 0, 0, 0, 0, 0, 0, t47, -t48, -t46, -g(1) * (t69 * t57 + t72) - g(2) * (t69 * t56 + t73) - g(3) * (t54 * pkin(4) - t55 * qJ(5) + t74);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:13
% EndTime: 2019-12-31 18:09:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->42), mult. (93->45), div. (0->0), fcn. (80->8), ass. (0->26)
t61 = qJ(2) + pkin(5);
t72 = g(3) * t61;
t64 = cos(qJ(3));
t50 = t64 * pkin(3) + pkin(2);
t59 = qJ(1) + pkin(7);
t52 = sin(t59);
t54 = cos(t59);
t63 = sin(qJ(1));
t56 = t63 * pkin(1);
t60 = -qJ(4) - pkin(6);
t71 = t52 * t50 + t54 * t60 + t56;
t62 = sin(qJ(3));
t70 = t62 * pkin(3) + t61;
t65 = cos(qJ(1));
t57 = t65 * pkin(1);
t69 = t54 * t50 - t52 * t60 + t57;
t68 = g(1) * t54 + g(2) * t52;
t67 = -g(1) * t65 - g(2) * t63;
t58 = qJ(3) + pkin(8);
t51 = sin(t58);
t53 = cos(t58);
t66 = pkin(4) * t53 + qJ(5) * t51;
t45 = g(1) * t52 - g(2) * t54;
t44 = -g(3) * t51 - t68 * t53;
t43 = -g(3) * t53 + t68 * t51;
t1 = [0, 0, 0, 0, 0, 0, t67, g(1) * t63 - g(2) * t65, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t68, t45, -g(3), t67 * pkin(1) - t72, 0, 0, 0, 0, 0, 0, -g(3) * t62 - t68 * t64, -g(3) * t64 + t68 * t62, -t45, -g(1) * (t54 * pkin(2) + t52 * pkin(6) + t57) - g(2) * (t52 * pkin(2) - t54 * pkin(6) + t56) - t72, 0, 0, 0, 0, 0, 0, t44, t43, -t45, -g(1) * t69 - g(2) * t71 - g(3) * t70, 0, 0, 0, 0, 0, 0, t44, -t45, -t43, -g(1) * (t66 * t54 + t69) - g(2) * (t66 * t52 + t71) - g(3) * (t51 * pkin(4) - t53 * qJ(5) + t70);];
U_reg = t1;

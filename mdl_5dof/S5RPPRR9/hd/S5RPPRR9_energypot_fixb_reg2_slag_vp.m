% Calculate inertial parameters regressor of potential energy for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:47
% EndTime: 2019-12-31 18:02:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (95->44), mult. (166->59), div. (0->0), fcn. (189->8), ass. (0->29)
t74 = g(3) * pkin(5);
t52 = -qJ(3) + pkin(5);
t73 = g(3) * t52;
t54 = sin(qJ(4));
t72 = g(3) * t54;
t71 = cos(qJ(1));
t70 = sin(qJ(1));
t53 = sin(qJ(5));
t56 = cos(qJ(4));
t69 = t53 * t56;
t55 = cos(qJ(5));
t68 = t55 * t56;
t67 = t71 * pkin(1) + t70 * qJ(2);
t66 = cos(pkin(8));
t65 = sin(pkin(8));
t64 = t71 * pkin(2) + t67;
t63 = -pkin(4) * t56 - pkin(7) * t54;
t40 = -t70 * t65 - t71 * t66;
t41 = t71 * t65 - t70 * t66;
t62 = g(1) * t41 - g(2) * t40;
t61 = g(1) * t40 + g(2) * t41;
t60 = t70 * pkin(1) - t71 * qJ(2);
t59 = t70 * pkin(2) + t60;
t58 = -t40 * pkin(3) + t41 * pkin(6) + t64;
t57 = -t41 * pkin(3) - t40 * pkin(6) + t59;
t43 = -g(1) * t71 - g(2) * t70;
t42 = g(1) * t70 - g(2) * t71;
t37 = g(3) * t56 - t61 * t54;
t1 = [0, 0, 0, 0, 0, 0, t43, t42, -g(3), -t74, 0, 0, 0, 0, 0, 0, t43, -g(3), -t42, -g(1) * t67 - g(2) * t60 - t74, 0, 0, 0, 0, 0, 0, t61, t62, g(3), -g(1) * t64 - g(2) * t59 - t73, 0, 0, 0, 0, 0, 0, t61 * t56 + t72, t37, -t62, -g(1) * t58 - g(2) * t57 - t73, 0, 0, 0, 0, 0, 0, -g(1) * (-t40 * t68 + t41 * t53) - g(2) * (-t40 * t53 - t41 * t68) + t55 * t72, -g(1) * (t40 * t69 + t41 * t55) - g(2) * (-t40 * t55 + t41 * t69) - t53 * t72, -t37, -g(1) * (t63 * t40 + t58) - g(2) * (t63 * t41 + t57) - g(3) * (-t54 * pkin(4) + t56 * pkin(7) + t52);];
U_reg = t1;

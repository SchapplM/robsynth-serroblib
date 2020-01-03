% Calculate inertial parameters regressor of potential energy for
% S5RPRRR13
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:26
% EndTime: 2019-12-31 19:15:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (87->58), mult. (132->71), div. (0->0), fcn. (127->8), ass. (0->35)
t57 = cos(qJ(4));
t45 = t57 * pkin(4) + pkin(3);
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t60 = -pkin(8) - pkin(7);
t80 = t45 * t55 + t58 * t60;
t79 = g(3) * pkin(5);
t78 = pkin(2) + pkin(5);
t59 = cos(qJ(1));
t77 = g(2) * t59;
t76 = g(3) * t58;
t53 = qJ(4) + qJ(5);
t46 = sin(t53);
t56 = sin(qJ(1));
t74 = t56 * t46;
t47 = cos(t53);
t73 = t56 * t47;
t54 = sin(qJ(4));
t72 = t56 * t54;
t71 = t56 * t57;
t69 = t59 * t46;
t68 = t59 * t47;
t67 = t59 * t54;
t66 = t59 * t57;
t49 = t56 * pkin(6);
t50 = t56 * pkin(1);
t65 = t49 + t50;
t64 = t59 * pkin(1) + t56 * qJ(2);
t63 = t59 * pkin(6) + t64;
t62 = -t59 * qJ(2) + t50;
t61 = pkin(3) * t55 - pkin(7) * t58;
t43 = g(1) * t56 - t77;
t44 = g(1) * t59 + g(2) * t56;
t42 = -g(3) * t55 + t43 * t58;
t1 = [0, 0, 0, 0, 0, 0, -t44, t43, -g(3), -t79, 0, 0, 0, 0, 0, 0, -g(3), t44, -t43, -g(1) * t64 - g(2) * t62 - t79, 0, 0, 0, 0, 0, 0, -t43 * t55 - t76, -t42, -t44, -g(1) * t63 - g(2) * (t49 + t62) - g(3) * t78, 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t71 + t67) - g(2) * (-t55 * t66 + t72) - t57 * t76, -g(1) * (-t55 * t72 + t66) - g(2) * (t55 * t67 + t71) + t54 * t76, t42, -g(1) * (t61 * t56 + t63) - g(2) * t65 - g(3) * (t58 * pkin(3) + t55 * pkin(7) + t78) - (-qJ(2) - t61) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t73 + t69) - g(2) * (-t55 * t68 + t74) - t47 * t76, -g(1) * (-t55 * t74 + t68) - g(2) * (t55 * t69 + t73) + t46 * t76, t42, -g(1) * (t80 * t56 + t63) - g(2) * (pkin(4) * t72 + t65) - g(3) * (t58 * t45 - t55 * t60 + t78) + (-g(1) * pkin(4) * t54 - g(2) * (-qJ(2) - t80)) * t59;];
U_reg = t1;

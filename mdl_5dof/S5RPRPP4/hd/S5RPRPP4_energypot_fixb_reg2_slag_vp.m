% Calculate inertial parameters regressor of potential energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:49
% EndTime: 2019-12-31 18:14:49
% DurationCPUTime: 0.07s
% Computational Cost: add. (81->43), mult. (99->44), div. (0->0), fcn. (86->6), ass. (0->24)
t51 = qJ(3) + pkin(7);
t45 = sin(t51);
t46 = cos(t51);
t68 = pkin(4) * t45 - qJ(5) * t46;
t67 = g(3) * pkin(5);
t66 = pkin(2) + pkin(5);
t53 = sin(qJ(3));
t65 = pkin(3) * t53;
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t63 = t56 * pkin(1) + t54 * qJ(2);
t55 = cos(qJ(3));
t61 = t55 * pkin(3) + t66;
t60 = t54 * t65 + t63;
t59 = -qJ(2) - t65;
t48 = t54 * pkin(1);
t52 = -qJ(4) - pkin(6);
t58 = -t54 * t52 + t48;
t57 = -t56 * qJ(2) + t48;
t42 = g(1) * t54 - g(2) * t56;
t43 = g(1) * t56 + g(2) * t54;
t41 = -g(3) * t45 + t42 * t46;
t40 = -g(3) * t46 - t42 * t45;
t1 = [0, 0, 0, 0, 0, 0, -t43, t42, -g(3), -t67, 0, 0, 0, 0, 0, 0, -g(3), t43, -t42, -g(1) * t63 - g(2) * t57 - t67, 0, 0, 0, 0, 0, 0, -g(3) * t55 - t42 * t53, g(3) * t53 - t42 * t55, -t43, -g(1) * (t56 * pkin(6) + t63) - g(2) * (t54 * pkin(6) + t57) - g(3) * t66, 0, 0, 0, 0, 0, 0, t40, -t41, -t43, -g(1) * (-t56 * t52 + t60) - g(2) * (t59 * t56 + t58) - g(3) * t61, 0, 0, 0, 0, 0, 0, t40, -t43, t41, -g(1) * (t68 * t54 + t60) - g(2) * t58 - g(3) * (t46 * pkin(4) + t45 * qJ(5) + t61) + (g(1) * t52 - g(2) * (t59 - t68)) * t56;];
U_reg = t1;

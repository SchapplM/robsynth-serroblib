% Calculate inertial parameters regressor of potential energy for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:29
% EndTime: 2019-12-31 18:47:29
% DurationCPUTime: 0.07s
% Computational Cost: add. (86->37), mult. (145->43), div. (0->0), fcn. (160->6), ass. (0->25)
t71 = g(3) * pkin(5);
t57 = -pkin(6) + pkin(5);
t70 = g(3) * t57;
t69 = cos(qJ(3));
t68 = sin(qJ(1));
t67 = sin(qJ(3));
t56 = cos(qJ(1));
t66 = t56 * pkin(1) + t68 * qJ(2);
t65 = t56 * pkin(2) + t66;
t64 = t68 * pkin(1) - t56 * qJ(2);
t63 = t68 * pkin(2) + t64;
t42 = -t56 * t69 - t68 * t67;
t43 = t56 * t67 - t68 * t69;
t62 = g(1) * t43 - g(2) * t42;
t61 = g(1) * t42 + g(2) * t43;
t54 = sin(qJ(4));
t55 = cos(qJ(4));
t60 = -pkin(4) * t55 - qJ(5) * t54;
t59 = -t42 * pkin(3) + t43 * pkin(7) + t65;
t58 = -t43 * pkin(3) - t42 * pkin(7) + t63;
t45 = -g(1) * t56 - g(2) * t68;
t44 = g(1) * t68 - g(2) * t56;
t38 = g(3) * t54 + t61 * t55;
t37 = g(3) * t55 - t61 * t54;
t1 = [0, 0, 0, 0, 0, 0, t45, t44, -g(3), -t71, 0, 0, 0, 0, 0, 0, t45, -g(3), -t44, -g(1) * t66 - g(2) * t64 - t71, 0, 0, 0, 0, 0, 0, t61, t62, g(3), -g(1) * t65 - g(2) * t63 - t70, 0, 0, 0, 0, 0, 0, t38, t37, -t62, -g(1) * t59 - g(2) * t58 - t70, 0, 0, 0, 0, 0, 0, t38, -t62, -t37, -g(1) * (t60 * t42 + t59) - g(2) * (t60 * t43 + t58) - g(3) * (-t54 * pkin(4) + t55 * qJ(5) + t57);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:00
% EndTime: 2020-01-03 11:34:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (118->41), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t66 = qJ(2) + pkin(5);
t54 = pkin(6) + t66;
t67 = g(1) * t54;
t56 = qJ(1) + pkin(8);
t49 = sin(t56);
t60 = sin(qJ(1));
t65 = t60 * pkin(1) + pkin(2) * t49;
t51 = cos(t56);
t61 = cos(qJ(1));
t64 = -t61 * pkin(1) - pkin(2) * t51;
t52 = qJ(3) + t56;
t45 = sin(t52);
t46 = cos(t52);
t63 = g(2) * t45 - g(3) * t46;
t62 = -g(2) * t60 + g(3) * t61;
t59 = -pkin(7) - qJ(4);
t58 = cos(pkin(9));
t57 = sin(pkin(9));
t55 = pkin(9) + qJ(5);
t50 = cos(t55);
t48 = sin(t55);
t47 = t58 * pkin(4) + pkin(3);
t43 = g(2) * t46 + g(3) * t45;
t1 = [0, 0, 0, 0, 0, 0, t62, -g(2) * t61 - g(3) * t60, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -g(2) * t49 + g(3) * t51, -g(2) * t51 - g(3) * t49, -g(1), t62 * pkin(1) - g(1) * t66, 0, 0, 0, 0, 0, 0, -t63, -t43, -g(1), -g(2) * t65 - g(3) * t64 - t67, 0, 0, 0, 0, 0, 0, -g(1) * t57 - t63 * t58, -g(1) * t58 + t63 * t57, t43, -t67 - g(2) * (t45 * pkin(3) - t46 * qJ(4) + t65) - g(3) * (-t46 * pkin(3) - t45 * qJ(4) + t64), 0, 0, 0, 0, 0, 0, -g(1) * t48 - t63 * t50, -g(1) * t50 + t63 * t48, t43, -g(1) * (t57 * pkin(4) + t54) - g(2) * (t45 * t47 + t46 * t59 + t65) - g(3) * (t45 * t59 - t46 * t47 + t64);];
U_reg = t1;

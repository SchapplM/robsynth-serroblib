% Calculate inertial parameters regressor of potential energy for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:28
% EndTime: 2020-01-03 11:45:28
% DurationCPUTime: 0.06s
% Computational Cost: add. (112->38), mult. (74->41), div. (0->0), fcn. (61->8), ass. (0->23)
t64 = qJ(2) + pkin(5);
t53 = pkin(6) + t64;
t65 = g(1) * t53;
t54 = qJ(1) + pkin(8);
t49 = sin(t54);
t57 = sin(qJ(1));
t63 = t57 * pkin(1) + pkin(2) * t49;
t50 = cos(t54);
t59 = cos(qJ(1));
t62 = -t59 * pkin(1) - pkin(2) * t50;
t51 = qJ(3) + t54;
t46 = sin(t51);
t47 = cos(t51);
t61 = g(2) * t46 - g(3) * t47;
t60 = -g(2) * t57 + g(3) * t59;
t58 = cos(qJ(4));
t56 = sin(qJ(4));
t55 = -qJ(5) - pkin(7);
t48 = t58 * pkin(4) + pkin(3);
t44 = g(2) * t47 + g(3) * t46;
t43 = -g(1) * t58 + t61 * t56;
t42 = -g(1) * t56 - t61 * t58;
t1 = [0, 0, 0, 0, 0, 0, t60, -g(2) * t59 - g(3) * t57, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -g(2) * t49 + g(3) * t50, -g(2) * t50 - g(3) * t49, -g(1), t60 * pkin(1) - g(1) * t64, 0, 0, 0, 0, 0, 0, -t61, -t44, -g(1), -g(2) * t63 - g(3) * t62 - t65, 0, 0, 0, 0, 0, 0, t42, t43, t44, -t65 - g(2) * (t46 * pkin(3) - t47 * pkin(7) + t63) - g(3) * (-t47 * pkin(3) - t46 * pkin(7) + t62), 0, 0, 0, 0, 0, 0, t42, t43, t44, -g(1) * (t56 * pkin(4) + t53) - g(2) * (t46 * t48 + t47 * t55 + t63) - g(3) * (t46 * t55 - t47 * t48 + t62);];
U_reg = t1;

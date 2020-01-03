% Calculate inertial parameters regressor of potential energy for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:42
% EndTime: 2019-12-31 16:30:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (63->38), mult. (128->50), div. (0->0), fcn. (130->6), ass. (0->26)
t48 = sin(qJ(2));
t50 = cos(qJ(2));
t64 = pkin(2) * t50 + pkin(5) * t48;
t61 = g(3) * t48;
t60 = g(3) * qJ(1);
t47 = sin(qJ(3));
t59 = t47 * t50;
t49 = cos(qJ(3));
t58 = t49 * t50;
t45 = sin(pkin(6));
t46 = cos(pkin(6));
t57 = t46 * pkin(1) + t45 * pkin(4);
t56 = t45 * pkin(1) - t46 * pkin(4);
t55 = t64 * t46 + t57;
t54 = t48 * pkin(2) - t50 * pkin(5) + qJ(1);
t53 = g(1) * t46 + g(2) * t45;
t52 = t64 * t45 + t56;
t30 = t45 * t59 + t46 * t49;
t32 = -t45 * t49 + t46 * t59;
t51 = g(1) * t32 + g(2) * t30 + t47 * t61;
t34 = g(1) * t45 - g(2) * t46;
t33 = t45 * t47 + t46 * t58;
t31 = t45 * t58 - t46 * t47;
t29 = -g(3) * t50 + t53 * t48;
t28 = -g(1) * t33 - g(2) * t31 - t49 * t61;
t1 = [0, 0, 0, 0, 0, 0, -t53, t34, -g(3), -t60, 0, 0, 0, 0, 0, 0, -t53 * t50 - t61, t29, -t34, -g(1) * t57 - g(2) * t56 - t60, 0, 0, 0, 0, 0, 0, t28, t51, -t29, -g(1) * t55 - g(2) * t52 - g(3) * t54, 0, 0, 0, 0, 0, 0, t28, -t29, -t51, -g(1) * (t33 * pkin(3) + t32 * qJ(4) + t55) - g(2) * (t31 * pkin(3) + t30 * qJ(4) + t52) - g(3) * ((pkin(3) * t49 + qJ(4) * t47) * t48 + t54);];
U_reg = t1;

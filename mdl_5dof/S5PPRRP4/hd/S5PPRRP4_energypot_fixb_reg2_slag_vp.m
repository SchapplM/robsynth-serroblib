% Calculate inertial parameters regressor of potential energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:40
% EndTime: 2019-12-31 17:34:40
% DurationCPUTime: 0.07s
% Computational Cost: add. (81->38), mult. (130->43), div. (0->0), fcn. (141->6), ass. (0->24)
t48 = -pkin(5) + qJ(1);
t61 = g(3) * t48;
t60 = cos(qJ(3));
t59 = sin(qJ(3));
t58 = g(3) * qJ(1);
t46 = cos(pkin(7));
t56 = sin(pkin(7));
t57 = t46 * pkin(1) + t56 * qJ(2);
t55 = t46 * pkin(2) + t57;
t54 = t56 * pkin(1) - t46 * qJ(2);
t53 = t56 * pkin(2) + t54;
t33 = -t46 * t60 - t56 * t59;
t34 = t46 * t59 - t56 * t60;
t52 = g(1) * t34 - g(2) * t33;
t51 = g(1) * t33 + g(2) * t34;
t50 = cos(qJ(4));
t49 = sin(qJ(4));
t47 = -qJ(5) - pkin(6);
t40 = t50 * pkin(4) + pkin(3);
t36 = -g(1) * t46 - g(2) * t56;
t35 = g(1) * t56 - g(2) * t46;
t31 = g(3) * t49 + t51 * t50;
t30 = g(3) * t50 - t51 * t49;
t1 = [0, 0, 0, 0, 0, 0, t36, t35, -g(3), -t58, 0, 0, 0, 0, 0, 0, t36, -g(3), -t35, -g(1) * t57 - g(2) * t54 - t58, 0, 0, 0, 0, 0, 0, t51, t52, g(3), -g(1) * t55 - g(2) * t53 - t61, 0, 0, 0, 0, 0, 0, t31, t30, -t52, -g(1) * (-t33 * pkin(3) + t34 * pkin(6) + t55) - g(2) * (-t34 * pkin(3) - t33 * pkin(6) + t53) - t61, 0, 0, 0, 0, 0, 0, t31, t30, -t52, -g(1) * (-t33 * t40 - t34 * t47 + t55) - g(2) * (t33 * t47 - t34 * t40 + t53) - g(3) * (-t49 * pkin(4) + t48);];
U_reg = t1;

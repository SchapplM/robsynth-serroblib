% Calculate inertial parameters regressor of potential energy for
% S5RPRRP1
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:56
% EndTime: 2019-12-05 17:59:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (77->43), mult. (92->42), div. (0->0), fcn. (79->6), ass. (0->22)
t65 = g(3) * pkin(5);
t64 = pkin(2) + pkin(5);
t59 = -pkin(7) - pkin(6);
t55 = sin(qJ(3));
t63 = t55 * pkin(3);
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t62 = t58 * pkin(1) + t56 * qJ(2);
t57 = cos(qJ(3));
t61 = t57 * pkin(3) + t64;
t50 = t56 * pkin(1);
t60 = -t58 * qJ(2) + t50;
t45 = g(1) * t56 - g(2) * t58;
t54 = qJ(3) + qJ(4);
t53 = -qJ(5) + t59;
t48 = cos(t54);
t47 = sin(t54);
t46 = g(1) * t58 + g(2) * t56;
t44 = pkin(4) * t47 + t63;
t43 = g(3) * t47 - t45 * t48;
t42 = -g(3) * t48 - t45 * t47;
t1 = [0, 0, 0, 0, 0, 0, -t46, t45, -g(3), -t65, 0, 0, 0, 0, 0, 0, -g(3), t46, -t45, -g(1) * t62 - g(2) * t60 - t65, 0, 0, 0, 0, 0, 0, -g(3) * t57 - t45 * t55, g(3) * t55 - t45 * t57, -t46, -g(1) * (t58 * pkin(6) + t62) - g(2) * (t56 * pkin(6) + t60) - g(3) * t64, 0, 0, 0, 0, 0, 0, t42, t43, -t46, -g(1) * (t56 * t63 - t58 * t59 + t62) - g(2) * (-t56 * t59 + t50 + (-qJ(2) - t63) * t58) - g(3) * t61, 0, 0, 0, 0, 0, 0, t42, t43, -t46, -g(1) * (t56 * t44 - t58 * t53 + t62) - g(2) * (-t56 * t53 + t50 + (-qJ(2) - t44) * t58) - g(3) * (pkin(4) * t48 + t61);];
U_reg = t1;

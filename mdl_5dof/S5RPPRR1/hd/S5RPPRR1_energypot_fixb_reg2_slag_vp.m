% Calculate inertial parameters regressor of potential energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:14
% EndTime: 2019-12-05 17:38:14
% DurationCPUTime: 0.06s
% Computational Cost: add. (62->40), mult. (80->38), div. (0->0), fcn. (67->6), ass. (0->21)
t54 = g(3) * pkin(5);
t53 = pkin(2) + pkin(5);
t42 = sin(qJ(4));
t52 = pkin(4) * t42;
t43 = sin(qJ(1));
t35 = t43 * qJ(3);
t38 = t43 * pkin(1);
t51 = t35 + t38;
t45 = cos(qJ(1));
t50 = t45 * pkin(1) + t43 * qJ(2);
t49 = pkin(3) + t53;
t48 = t45 * qJ(3) + t50;
t47 = -t45 * qJ(2) + t38;
t32 = g(1) * t45 + g(2) * t43;
t46 = -pkin(7) - pkin(6);
t44 = cos(qJ(4));
t41 = qJ(4) + qJ(5);
t34 = cos(t41);
t33 = sin(t41);
t31 = g(1) * t43 - g(2) * t45;
t1 = [0, 0, 0, 0, 0, 0, -t32, t31, -g(3), -t54, 0, 0, 0, 0, 0, 0, -g(3), t32, -t31, -g(1) * t50 - g(2) * t47 - t54, 0, 0, 0, 0, 0, 0, -g(3), -t31, -t32, -g(1) * t48 - g(2) * (t35 + t47) - g(3) * t53, 0, 0, 0, 0, 0, 0, -g(3) * t44 - t32 * t42, g(3) * t42 - t32 * t44, t31, -g(1) * (-t43 * pkin(6) + t48) - g(2) * ((pkin(6) - qJ(2)) * t45 + t51) - g(3) * t49, 0, 0, 0, 0, 0, 0, -g(3) * t34 - t32 * t33, g(3) * t33 - t32 * t34, t31, -g(1) * (t43 * t46 + t45 * t52 + t48) - g(2) * (t43 * t52 + (-qJ(2) - t46) * t45 + t51) - g(3) * (t44 * pkin(4) + t49);];
U_reg = t1;

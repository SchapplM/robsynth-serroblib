% Calculate inertial parameters regressor of potential energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:21
% EndTime: 2019-12-05 15:03:21
% DurationCPUTime: 0.09s
% Computational Cost: add. (113->52), mult. (130->62), div. (0->0), fcn. (121->8), ass. (0->32)
t45 = cos(pkin(7));
t41 = pkin(8) + qJ(3);
t37 = sin(t41);
t54 = qJ(4) * t37;
t38 = cos(t41);
t60 = t38 * t45;
t64 = pkin(3) * t60 + t45 * t54;
t63 = g(3) * t38;
t62 = g(3) * qJ(1);
t43 = sin(pkin(7));
t61 = t38 * t43;
t47 = sin(qJ(5));
t59 = t43 * t47;
t48 = cos(qJ(5));
t58 = t43 * t48;
t57 = t45 * t47;
t56 = t45 * t48;
t44 = cos(pkin(8));
t36 = t44 * pkin(2) + pkin(1);
t46 = -pkin(5) - qJ(2);
t55 = t43 * t36 + t45 * t46;
t42 = sin(pkin(8));
t53 = t42 * pkin(2) + qJ(1);
t31 = t45 * t36;
t52 = -t43 * t46 + t31;
t51 = pkin(3) * t61 + t43 * t54 + t55;
t50 = g(1) * t45 + g(2) * t43;
t49 = t37 * pkin(3) - t38 * qJ(4) + t53;
t27 = g(1) * t43 - g(2) * t45;
t26 = g(3) * t37 + t50 * t38;
t25 = t50 * t37 - t63;
t1 = [0, 0, 0, 0, 0, 0, -t50, t27, -g(3), -t62, 0, 0, 0, 0, 0, 0, -g(3) * t42 - t50 * t44, -g(3) * t44 + t50 * t42, -t27, -g(1) * (t45 * pkin(1) + t43 * qJ(2)) - g(2) * (t43 * pkin(1) - t45 * qJ(2)) - t62, 0, 0, 0, 0, 0, 0, -t26, t25, -t27, -g(1) * t52 - g(2) * t55 - g(3) * t53, 0, 0, 0, 0, 0, 0, -t27, t26, -t25, -g(1) * (t52 + t64) - g(2) * t51 - g(3) * t49, 0, 0, 0, 0, 0, 0, -g(1) * (t37 * t57 + t58) - g(2) * (t37 * t59 - t56) + t47 * t63, -g(1) * (t37 * t56 - t59) - g(2) * (t37 * t58 + t57) + t48 * t63, -t26, -g(1) * (pkin(6) * t60 + t31 + (pkin(4) - t46) * t43 + t64) - g(2) * (-t45 * pkin(4) + pkin(6) * t61 + t51) - g(3) * (t37 * pkin(6) + t49);];
U_reg = t1;

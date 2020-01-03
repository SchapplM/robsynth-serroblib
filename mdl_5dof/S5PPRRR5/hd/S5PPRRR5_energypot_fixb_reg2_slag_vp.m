% Calculate inertial parameters regressor of potential energy for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:44
% EndTime: 2019-12-31 17:35:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (101->42), mult. (115->51), div. (0->0), fcn. (120->8), ass. (0->28)
t57 = -pkin(5) + qJ(1);
t60 = g(3) * (-pkin(6) + t57);
t59 = g(3) * qJ(1);
t42 = sin(pkin(8));
t45 = sin(qJ(3));
t58 = t42 * t45;
t43 = cos(pkin(8));
t56 = t43 * pkin(1) + t42 * qJ(2);
t55 = qJ(3) + qJ(4);
t54 = cos(t55);
t53 = sin(t55);
t39 = t42 * pkin(1);
t52 = -t43 * qJ(2) + t39;
t47 = cos(qJ(3));
t37 = t47 * pkin(3) + pkin(2);
t51 = pkin(3) * t58 + t43 * t37 + t56;
t25 = -t42 * t53 - t43 * t54;
t26 = -t42 * t54 + t43 * t53;
t50 = g(1) * t26 - g(2) * t25;
t49 = g(1) * t25 + g(2) * t26;
t48 = t42 * t37 + t39 + (-pkin(3) * t45 - qJ(2)) * t43;
t46 = cos(qJ(5));
t44 = sin(qJ(5));
t30 = -g(1) * t43 - g(2) * t42;
t29 = g(1) * t42 - g(2) * t43;
t28 = t42 * t47 - t43 * t45;
t27 = -t43 * t47 - t58;
t1 = [0, 0, 0, 0, 0, 0, t30, t29, -g(3), -t59, 0, 0, 0, 0, 0, 0, t30, -g(3), -t29, -g(1) * t56 - g(2) * t52 - t59, 0, 0, 0, 0, 0, 0, g(1) * t27 - g(2) * t28, -g(1) * t28 - g(2) * t27, g(3), -g(1) * (t43 * pkin(2) + t56) - g(2) * (t42 * pkin(2) + t52) - g(3) * t57, 0, 0, 0, 0, 0, 0, t49, t50, g(3), -g(1) * t51 - g(2) * t48 - t60, 0, 0, 0, 0, 0, 0, g(3) * t44 + t49 * t46, g(3) * t46 - t49 * t44, -t50, -g(1) * (-t25 * pkin(4) + t26 * pkin(7) + t51) - g(2) * (-t26 * pkin(4) - t25 * pkin(7) + t48) - t60;];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:41
% EndTime: 2019-12-31 18:16:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (64->39), mult. (106->42), div. (0->0), fcn. (93->4), ass. (0->22)
t60 = g(3) * pkin(5);
t59 = pkin(2) + pkin(5);
t46 = sin(qJ(3));
t47 = sin(qJ(1));
t58 = t46 * t47;
t49 = cos(qJ(1));
t57 = t49 * pkin(1) + t47 * qJ(2);
t48 = cos(qJ(3));
t56 = qJ(4) * t48;
t41 = t47 * pkin(6);
t42 = t47 * pkin(1);
t55 = t49 * t56 + t41 + t42;
t54 = t49 * pkin(6) + t57;
t53 = t48 * pkin(3) + t46 * qJ(4) + t59;
t52 = -pkin(3) * t46 - qJ(2);
t51 = -t49 * qJ(2) + t42;
t35 = g(1) * t47 - g(2) * t49;
t50 = pkin(3) * t58 - t47 * t56 + t54;
t36 = g(1) * t49 + g(2) * t47;
t34 = -g(3) * t46 + t35 * t48;
t33 = -g(3) * t48 - t35 * t46;
t1 = [0, 0, 0, 0, 0, 0, -t36, t35, -g(3), -t60, 0, 0, 0, 0, 0, 0, -g(3), t36, -t35, -g(1) * t57 - g(2) * t51 - t60, 0, 0, 0, 0, 0, 0, t33, -t34, -t36, -g(1) * t54 - g(2) * (t41 + t51) - g(3) * t59, 0, 0, 0, 0, 0, 0, t33, -t36, t34, -g(1) * t50 - g(2) * (t52 * t49 + t55) - g(3) * t53, 0, 0, 0, 0, 0, 0, t33, t34, t36, -g(1) * (pkin(4) * t58 + t50) - g(2) * (-t47 * qJ(5) + t55) - g(3) * (t48 * pkin(4) + t53) + (g(1) * qJ(5) - g(2) * (-pkin(4) * t46 + t52)) * t49;];
U_reg = t1;

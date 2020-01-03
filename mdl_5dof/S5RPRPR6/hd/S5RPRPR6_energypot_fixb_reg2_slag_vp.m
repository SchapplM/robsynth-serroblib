% Calculate inertial parameters regressor of potential energy for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:48
% EndTime: 2019-12-31 18:17:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (108->36), mult. (65->37), div. (0->0), fcn. (52->8), ass. (0->21)
t55 = qJ(2) + pkin(5);
t44 = pkin(6) + t55;
t56 = g(3) * t44;
t45 = qJ(1) + pkin(8);
t39 = sin(t45);
t47 = sin(qJ(1));
t54 = t47 * pkin(1) + pkin(2) * t39;
t40 = cos(t45);
t49 = cos(qJ(1));
t53 = t49 * pkin(1) + pkin(2) * t40;
t41 = qJ(3) + t45;
t37 = sin(t41);
t38 = cos(t41);
t52 = t38 * pkin(3) + t37 * qJ(4) + t53;
t30 = g(1) * t37 - g(2) * t38;
t51 = -g(1) * t49 - g(2) * t47;
t50 = t37 * pkin(3) - t38 * qJ(4) + t54;
t48 = cos(qJ(5));
t46 = sin(qJ(5));
t31 = g(1) * t38 + g(2) * t37;
t1 = [0, 0, 0, 0, 0, 0, t51, g(1) * t47 - g(2) * t49, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t40 - g(2) * t39, g(1) * t39 - g(2) * t40, -g(3), t51 * pkin(1) - g(3) * t55, 0, 0, 0, 0, 0, 0, -t31, t30, -g(3), -g(1) * t53 - g(2) * t54 - t56, 0, 0, 0, 0, 0, 0, -g(3), t31, -t30, -g(1) * t52 - g(2) * t50 - t56, 0, 0, 0, 0, 0, 0, -g(3) * t48 - t30 * t46, g(3) * t46 - t30 * t48, -t31, -g(1) * (t38 * pkin(7) + t52) - g(2) * (t37 * pkin(7) + t50) - g(3) * (pkin(4) + t44);];
U_reg = t1;

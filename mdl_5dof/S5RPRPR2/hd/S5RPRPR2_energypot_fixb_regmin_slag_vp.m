% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:06
% EndTime: 2022-01-23 09:19:06
% DurationCPUTime: 0.04s
% Computational Cost: add. (64->23), mult. (41->27), div. (0->0), fcn. (35->10), ass. (0->14)
t54 = qJ(2) + pkin(5);
t48 = qJ(1) + pkin(8);
t46 = qJ(3) + t48;
t42 = sin(t46);
t43 = cos(t46);
t53 = g(1) * t43 + g(2) * t42;
t50 = sin(qJ(1));
t51 = cos(qJ(1));
t52 = -g(1) * t51 - g(2) * t50;
t47 = pkin(9) + qJ(5);
t45 = cos(t47);
t44 = sin(t47);
t41 = g(1) * t42 - g(2) * t43;
t1 = [0, t52, g(1) * t50 - g(2) * t51, t52 * pkin(1) - g(3) * t54, 0, -t53, t41, -g(3) * sin(pkin(9)) - t53 * cos(pkin(9)), -t41, -g(1) * (t43 * pkin(3) + t42 * qJ(4) + pkin(2) * cos(t48) + t51 * pkin(1)) - g(2) * (t42 * pkin(3) - t43 * qJ(4) + pkin(2) * sin(t48) + t50 * pkin(1)) - g(3) * (pkin(6) + t54), 0, 0, 0, 0, 0, -g(3) * t44 - t53 * t45, -g(3) * t45 + t53 * t44;];
U_reg = t1;

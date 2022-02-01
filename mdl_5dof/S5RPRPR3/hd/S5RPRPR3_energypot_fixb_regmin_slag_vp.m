% Calculate minimal parameter regressor of potential energy for
% S5RPRPR3
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:54
% EndTime: 2022-01-23 09:20:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (70->27), mult. (51->39), div. (0->0), fcn. (49->10), ass. (0->17)
t61 = g(3) * sin(pkin(9));
t51 = cos(pkin(9));
t52 = sin(qJ(5));
t60 = t51 * t52;
t54 = cos(qJ(5));
t59 = t51 * t54;
t58 = qJ(2) + pkin(5);
t49 = qJ(1) + pkin(8);
t48 = qJ(3) + t49;
t46 = sin(t48);
t47 = cos(t48);
t57 = -g(1) * t47 - g(2) * t46;
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t56 = -g(1) * t55 - g(2) * t53;
t45 = g(1) * t46 - g(2) * t47;
t1 = [0, t56, g(1) * t53 - g(2) * t55, t56 * pkin(1) - g(3) * t58, 0, t57, t45, t57 * t51 - t61, -t45, -g(1) * (t47 * pkin(3) + t46 * qJ(4) + pkin(2) * cos(t49) + t55 * pkin(1)) - g(2) * (t46 * pkin(3) - t47 * qJ(4) + pkin(2) * sin(t49) + t53 * pkin(1)) - g(3) * (pkin(6) + t58), 0, 0, 0, 0, 0, -g(1) * (t46 * t52 + t47 * t59) - g(2) * (t46 * t59 - t47 * t52) - t54 * t61, -g(1) * (t46 * t54 - t47 * t60) - g(2) * (-t46 * t60 - t47 * t54) + t52 * t61;];
U_reg = t1;

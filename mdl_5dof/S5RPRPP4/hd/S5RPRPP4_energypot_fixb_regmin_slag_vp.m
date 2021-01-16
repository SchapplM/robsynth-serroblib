% Calculate minimal parameter regressor of potential energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:00
% EndTime: 2021-01-15 11:24:00
% DurationCPUTime: 0.05s
% Computational Cost: add. (70->35), mult. (81->41), div. (0->0), fcn. (72->8), ass. (0->20)
t59 = cos(qJ(3));
t62 = t59 * pkin(3) + pkin(2) + pkin(5);
t58 = sin(qJ(1));
t60 = cos(qJ(1));
t45 = g(1) * t58 - g(2) * t60;
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t57 = sin(qJ(3));
t61 = (pkin(4) * t56 + qJ(5) * t55 + pkin(3)) * t57 - (-t55 * pkin(4) + qJ(5) * t56) * t59 + qJ(2);
t54 = qJ(3) + pkin(7);
t53 = pkin(1) + pkin(6) + qJ(4);
t51 = cos(t54);
t50 = sin(t54);
t49 = t57 * pkin(3) + qJ(2);
t48 = t53 * t60;
t47 = t53 * t58;
t46 = g(1) * t60 + g(2) * t58;
t42 = -g(3) * t50 + t45 * t51;
t41 = -g(3) * t51 - t45 * t50;
t1 = [0, -t46, t45, t46, -t45, -g(1) * (t60 * pkin(1) + t58 * qJ(2)) - g(2) * (t58 * pkin(1) - t60 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t59 - t45 * t57, g(3) * t57 - t45 * t59, t41, -t42, -t46, -g(1) * (t49 * t58 + t48) - g(2) * (-t49 * t60 + t47) - g(3) * t62, t41, -t46, t42, -g(1) * (t61 * t58 + t48) - g(2) * (-t61 * t60 + t47) - g(3) * (t51 * pkin(4) + t50 * qJ(5) + t62);];
U_reg = t1;

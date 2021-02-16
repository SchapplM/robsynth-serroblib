% Calculate minimal parameter regressor of potential energy for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:51
% EndTime: 2021-01-15 11:33:51
% DurationCPUTime: 0.04s
% Computational Cost: add. (51->25), mult. (57->32), div. (0->0), fcn. (51->8), ass. (0->15)
t56 = qJ(3) + pkin(8);
t58 = sin(qJ(1));
t60 = cos(qJ(1));
t47 = g(1) * t58 - g(2) * t60;
t59 = cos(qJ(3));
t57 = sin(qJ(3));
t55 = pkin(1) + pkin(6) + qJ(4);
t54 = qJ(5) + t56;
t53 = cos(t56);
t52 = sin(t56);
t51 = t57 * pkin(3) + qJ(2);
t50 = cos(t54);
t49 = sin(t54);
t48 = g(1) * t60 + g(2) * t58;
t1 = [0, -t48, t47, t48, -t47, -g(1) * (t60 * pkin(1) + t58 * qJ(2)) - g(2) * (t58 * pkin(1) - t60 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t59 - t47 * t57, g(3) * t57 - t47 * t59, -g(3) * t53 - t47 * t52, g(3) * t52 - t47 * t53, -t48, -g(1) * (t51 * t58 + t55 * t60) - g(2) * (-t51 * t60 + t55 * t58) - g(3) * (t59 * pkin(3) + pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t50 - t47 * t49, g(3) * t49 - t47 * t50;];
U_reg = t1;

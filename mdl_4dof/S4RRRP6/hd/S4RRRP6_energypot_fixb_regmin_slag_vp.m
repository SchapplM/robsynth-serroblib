% Calculate minimal parameter regressor of potential energy for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:13
% EndTime: 2021-01-15 14:39:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (42->25), mult. (77->36), div. (0->0), fcn. (82->6), ass. (0->19)
t48 = sin(qJ(2));
t57 = g(3) * t48;
t49 = sin(qJ(1));
t51 = cos(qJ(2));
t56 = t49 * t51;
t47 = sin(qJ(3));
t52 = cos(qJ(1));
t55 = t52 * t47;
t50 = cos(qJ(3));
t54 = t52 * t50;
t53 = g(1) * t52 + g(2) * t49;
t46 = qJ(4) + pkin(6);
t45 = t50 * pkin(3) + pkin(2);
t44 = t47 * pkin(3) + pkin(5);
t43 = t45 * t51 + t46 * t48 + pkin(1);
t42 = -g(3) * t51 + t53 * t48;
t41 = -g(1) * (t49 * t47 + t51 * t54) - g(2) * (t50 * t56 - t55) - t50 * t57;
t40 = -g(1) * (t49 * t50 - t51 * t55) - g(2) * (-t47 * t56 - t54) + t47 * t57;
t1 = [0, -t53, g(1) * t49 - g(2) * t52, 0, 0, 0, 0, 0, -t53 * t51 - t57, t42, 0, 0, 0, 0, 0, t41, t40, t41, t40, -t42, -g(1) * (t43 * t52 + t44 * t49) - g(2) * (t43 * t49 - t44 * t52) - g(3) * (t48 * t45 - t51 * t46 + pkin(4));];
U_reg = t1;

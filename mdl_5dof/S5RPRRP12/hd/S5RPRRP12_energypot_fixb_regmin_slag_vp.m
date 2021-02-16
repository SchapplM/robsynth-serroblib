% Calculate minimal parameter regressor of potential energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:52
% EndTime: 2021-01-15 19:25:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (51->32), mult. (88->42), div. (0->0), fcn. (90->6), ass. (0->21)
t48 = cos(qJ(3));
t55 = g(3) * t48;
t44 = sin(qJ(4));
t46 = sin(qJ(1));
t54 = t46 * t44;
t47 = cos(qJ(4));
t53 = t46 * t47;
t49 = cos(qJ(1));
t52 = t49 * t44;
t51 = t49 * t47;
t39 = g(1) * t46 - g(2) * t49;
t42 = t47 * pkin(4) + pkin(3);
t43 = -qJ(5) - pkin(7);
t45 = sin(qJ(3));
t50 = t42 * t45 + t43 * t48 + qJ(2);
t41 = t44 * pkin(4) + pkin(1) + pkin(6);
t40 = g(1) * t49 + g(2) * t46;
t38 = -g(3) * t45 + t39 * t48;
t37 = -g(1) * (t45 * t53 + t52) - g(2) * (-t45 * t51 + t54) - t47 * t55;
t36 = -g(1) * (-t45 * t54 + t51) - g(2) * (t45 * t52 + t53) + t44 * t55;
t1 = [0, -t40, t39, t40, -t39, -g(1) * (t49 * pkin(1) + t46 * qJ(2)) - g(2) * (t46 * pkin(1) - t49 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t39 * t45 - t55, -t38, 0, 0, 0, 0, 0, t37, t36, t37, t36, t38, -g(1) * (t41 * t49 + t50 * t46) - g(2) * (t41 * t46 - t50 * t49) - g(3) * (t48 * t42 - t43 * t45 + pkin(2) + pkin(5));];
U_reg = t1;

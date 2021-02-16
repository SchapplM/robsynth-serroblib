% Calculate minimal parameter regressor of potential energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:08
% EndTime: 2021-01-15 15:42:08
% DurationCPUTime: 0.04s
% Computational Cost: add. (68->24), mult. (49->28), div. (0->0), fcn. (45->10), ass. (0->16)
t57 = qJ(3) + pkin(9);
t56 = pkin(8) + qJ(2);
t51 = sin(t56);
t53 = cos(t56);
t61 = g(1) * t53 + g(2) * t51;
t60 = cos(qJ(3));
t59 = sin(qJ(3));
t58 = -pkin(6) - qJ(4);
t55 = qJ(5) + t57;
t54 = cos(t57);
t52 = sin(t57);
t50 = t60 * pkin(3) + pkin(2);
t49 = cos(t55);
t48 = sin(t55);
t47 = g(1) * t51 - g(2) * t53;
t1 = [-g(3) * qJ(1), 0, -t61, t47, 0, 0, 0, 0, 0, -g(3) * t59 - t61 * t60, -g(3) * t60 + t61 * t59, -g(3) * t52 - t61 * t54, -g(3) * t54 + t61 * t52, -t47, -g(1) * (t53 * t50 - t51 * t58 + cos(pkin(8)) * pkin(1)) - g(2) * (t51 * t50 + t53 * t58 + sin(pkin(8)) * pkin(1)) - g(3) * (t59 * pkin(3) + pkin(5) + qJ(1)), 0, 0, 0, 0, 0, -g(3) * t48 - t61 * t49, -g(3) * t49 + t61 * t48;];
U_reg = t1;

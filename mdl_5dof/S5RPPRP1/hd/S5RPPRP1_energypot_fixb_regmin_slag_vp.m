% Calculate minimal parameter regressor of potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:55:46
% EndTime: 2021-01-15 16:55:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (90->36), mult. (96->53), div. (0->0), fcn. (95->8), ass. (0->25)
t53 = sin(pkin(8));
t70 = g(1) * t53;
t56 = qJ(2) + pkin(5);
t69 = g(1) * t56;
t60 = cos(qJ(1));
t68 = t60 * pkin(1);
t54 = cos(pkin(8));
t57 = sin(qJ(4));
t67 = t54 * t57;
t59 = cos(qJ(4));
t66 = t54 * t59;
t52 = qJ(1) + pkin(7);
t49 = sin(t52);
t58 = sin(qJ(1));
t65 = t58 * pkin(1) + t49 * pkin(2);
t64 = pkin(4) * t57 + qJ(3);
t50 = cos(t52);
t63 = -g(2) * t49 + g(3) * t50;
t62 = -g(2) * t58 + g(3) * t60;
t48 = t59 * pkin(4) + pkin(3);
t55 = -qJ(5) - pkin(6);
t61 = t48 * t54 - t53 * t55;
t46 = -t59 * t70 - g(2) * (t49 * t66 - t50 * t57) - g(3) * (-t49 * t57 - t50 * t66);
t45 = t57 * t70 - g(2) * (-t49 * t67 - t50 * t59) - g(3) * (-t49 * t59 + t50 * t67);
t1 = [0, t62, -g(2) * t60 - g(3) * t58, t62 * pkin(1) - t69, t63 * t54 - t70, g(2) * t50 + g(3) * t49, -t69 - g(2) * (-t50 * qJ(3) + t65) - g(3) * (-t50 * pkin(2) - t49 * qJ(3) - t68), 0, 0, 0, 0, 0, t46, t45, t46, t45, g(1) * t54 + t63 * t53, -g(1) * (t53 * t48 + t54 * t55 + t56) - g(2) * t65 + g(3) * t68 + (-g(2) * t61 + g(3) * t64) * t49 + (g(2) * t64 - g(3) * (-pkin(2) - t61)) * t50;];
U_reg = t1;

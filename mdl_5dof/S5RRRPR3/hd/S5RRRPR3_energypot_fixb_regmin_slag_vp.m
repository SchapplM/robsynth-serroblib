% Calculate minimal parameter regressor of potential energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:59:34
% EndTime: 2021-01-15 22:59:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (69->25), mult. (52->31), div. (0->0), fcn. (49->10), ass. (0->18)
t62 = qJ(3) + pkin(9);
t63 = qJ(1) + qJ(2);
t60 = sin(t63);
t61 = cos(t63);
t69 = g(2) * t60 - g(3) * t61;
t68 = cos(qJ(1));
t67 = cos(qJ(3));
t66 = sin(qJ(1));
t65 = sin(qJ(3));
t64 = -qJ(4) - pkin(7);
t59 = qJ(5) + t62;
t58 = cos(t62);
t57 = sin(t62);
t56 = t67 * pkin(3) + pkin(2);
t55 = cos(t59);
t54 = sin(t59);
t53 = g(2) * t61 + g(3) * t60;
t1 = [0, -g(2) * t66 + g(3) * t68, -g(2) * t68 - g(3) * t66, 0, -t69, -t53, 0, 0, 0, 0, 0, -g(1) * t65 - t69 * t67, -g(1) * t67 + t69 * t65, -g(1) * t57 - t69 * t58, -g(1) * t58 + t69 * t57, t53, -g(1) * (t65 * pkin(3) + pkin(5) + pkin(6)) - g(2) * (t66 * pkin(1) + t60 * t56 + t61 * t64) - g(3) * (-t68 * pkin(1) - t61 * t56 + t60 * t64), 0, 0, 0, 0, 0, -g(1) * t54 - t69 * t55, -g(1) * t55 + t69 * t54;];
U_reg = t1;

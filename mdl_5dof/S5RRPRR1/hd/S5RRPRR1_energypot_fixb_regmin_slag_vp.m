% Calculate minimal parameter regressor of potential energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:17
% EndTime: 2021-01-15 21:13:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (39->21), mult. (66->34), div. (0->0), fcn. (67->8), ass. (0->21)
t56 = cos(qJ(2));
t65 = pkin(1) * t56;
t51 = qJ(2) + qJ(4);
t49 = sin(t51);
t64 = g(3) * t49;
t53 = sin(qJ(2));
t63 = g(3) * t53;
t52 = sin(qJ(5));
t54 = sin(qJ(1));
t62 = t54 * t52;
t55 = cos(qJ(5));
t61 = t54 * t55;
t57 = cos(qJ(1));
t60 = t57 * t52;
t59 = t57 * t55;
t58 = g(1) * t57 + g(2) * t54;
t50 = cos(t51);
t48 = g(1) * t54 - g(2) * t57;
t47 = -t58 * t56 - t63;
t46 = -g(3) * t56 + t58 * t53;
t1 = [0, -t58, t48, 0, 0, 0, 0, 0, t47, t46, t47, t46, -t48, -g(1) * (t54 * qJ(3) + t57 * t65) - g(2) * (-t57 * qJ(3) + t54 * t65) - pkin(1) * t63, 0, 0, 0, 0, 0, -t58 * t50 - t64, -g(3) * t50 + t58 * t49, 0, 0, 0, 0, 0, -g(1) * (t50 * t59 + t62) - g(2) * (t50 * t61 - t60) - t55 * t64, -g(1) * (-t50 * t60 + t61) - g(2) * (-t50 * t62 - t59) + t52 * t64;];
U_reg = t1;

% Calculate minimal parameter regressor of potential energy for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:46
% EndTime: 2021-01-15 22:14:46
% DurationCPUTime: 0.05s
% Computational Cost: add. (98->31), mult. (76->36), div. (0->0), fcn. (70->8), ass. (0->21)
t63 = sin(qJ(3));
t71 = t63 * pkin(3) + pkin(5) + pkin(6);
t65 = cos(qJ(3));
t52 = t65 * pkin(3) + pkin(2);
t61 = qJ(1) + qJ(2);
t55 = sin(t61);
t56 = cos(t61);
t62 = -qJ(4) - pkin(7);
t64 = sin(qJ(1));
t70 = t64 * pkin(1) + t55 * t52 + t56 * t62;
t66 = cos(qJ(1));
t69 = t66 * pkin(1) + t56 * t52 - t55 * t62;
t68 = g(1) * t56 + g(2) * t55;
t60 = qJ(3) + pkin(8);
t53 = sin(t60);
t54 = cos(t60);
t67 = pkin(4) * t54 + qJ(5) * t53;
t47 = g(1) * t55 - g(2) * t56;
t46 = -g(3) * t53 - t68 * t54;
t45 = -g(3) * t54 + t68 * t53;
t1 = [0, -g(1) * t66 - g(2) * t64, g(1) * t64 - g(2) * t66, 0, -t68, t47, 0, 0, 0, 0, 0, -g(3) * t63 - t68 * t65, -g(3) * t65 + t68 * t63, t46, t45, -t47, -g(1) * t69 - g(2) * t70 - g(3) * t71, t46, -t47, -t45, -g(1) * (t67 * t56 + t69) - g(2) * (t67 * t55 + t70) - g(3) * (t53 * pkin(4) - t54 * qJ(5) + t71);];
U_reg = t1;

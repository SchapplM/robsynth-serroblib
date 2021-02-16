% Calculate minimal parameter regressor of potential energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:18
% EndTime: 2021-01-15 20:52:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (65->31), mult. (104->40), div. (0->0), fcn. (110->6), ass. (0->20)
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t64 = g(1) * t63 + g(2) * t60;
t62 = cos(qJ(2));
t61 = cos(qJ(4));
t59 = sin(qJ(2));
t58 = sin(qJ(4));
t57 = qJ(5) - pkin(6) + pkin(7);
t54 = t58 * pkin(4) + qJ(3);
t53 = t61 * pkin(4) + pkin(2) + pkin(3);
t52 = g(1) * t60 - g(2) * t63;
t51 = t62 * pkin(2) + t59 * qJ(3) + pkin(1);
t50 = -t62 * t58 + t59 * t61;
t49 = -t59 * t58 - t62 * t61;
t48 = -g(3) * t59 - t64 * t62;
t47 = -g(3) * t62 + t64 * t59;
t46 = t53 * t62 + t54 * t59 + pkin(1);
t45 = -g(3) * t50 + t64 * t49;
t44 = -g(3) * t49 - t64 * t50;
t1 = [0, -t64, t52, 0, 0, 0, 0, 0, t48, t47, t48, -t52, -t47, -g(1) * (t60 * pkin(6) + t51 * t63) - g(2) * (-t63 * pkin(6) + t51 * t60) - g(3) * (t59 * pkin(2) - t62 * qJ(3) + pkin(5)), 0, 0, 0, 0, 0, t45, t44, t45, t44, t52, -g(1) * (t46 * t63 - t57 * t60) - g(2) * (t46 * t60 + t57 * t63) - g(3) * (t53 * t59 - t54 * t62 + pkin(5));];
U_reg = t1;

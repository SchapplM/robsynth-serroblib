% Calculate minimal parameter regressor of potential energy for
% S5RRPRP10
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
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:44
% EndTime: 2021-01-15 21:02:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (62->37), mult. (102->46), div. (0->0), fcn. (104->6), ass. (0->23)
t58 = cos(qJ(2));
t65 = g(3) * t58;
t54 = sin(qJ(4));
t56 = sin(qJ(1));
t64 = t56 * t54;
t57 = cos(qJ(4));
t63 = t56 * t57;
t59 = cos(qJ(1));
t62 = t59 * t54;
t61 = t59 * t57;
t60 = g(1) * t59 + g(2) * t56;
t55 = sin(qJ(2));
t53 = pkin(2) + pkin(7) + qJ(5);
t52 = t54 * pkin(4) + qJ(3);
t51 = t57 * pkin(4) + pkin(3) + pkin(6);
t50 = g(1) * t56 - g(2) * t59;
t49 = t58 * pkin(2) + t55 * qJ(3) + pkin(1);
t48 = t52 * t55 + t53 * t58 + pkin(1);
t47 = g(3) * t55 + t60 * t58;
t46 = t60 * t55 - t65;
t45 = -g(1) * (t55 * t62 + t63) - g(2) * (t55 * t64 - t61) + t54 * t65;
t44 = -g(1) * (t55 * t61 - t64) - g(2) * (t55 * t63 + t62) + t57 * t65;
t1 = [0, -t60, t50, 0, 0, 0, 0, 0, -t47, t46, -t50, t47, -t46, -g(1) * (t56 * pkin(6) + t49 * t59) - g(2) * (-t59 * pkin(6) + t49 * t56) - g(3) * (t55 * pkin(2) - t58 * qJ(3) + pkin(5)), 0, 0, 0, 0, 0, t45, t44, t45, t44, -t47, -g(1) * (t48 * t59 + t51 * t56) - g(2) * (t48 * t56 - t51 * t59) - g(3) * (-t52 * t58 + t53 * t55 + pkin(5));];
U_reg = t1;

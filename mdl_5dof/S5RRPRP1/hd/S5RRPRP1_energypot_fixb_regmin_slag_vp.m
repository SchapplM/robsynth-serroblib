% Calculate minimal parameter regressor of potential energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:51
% EndTime: 2021-01-15 20:08:51
% DurationCPUTime: 0.11s
% Computational Cost: add. (70->24), mult. (51->32), div. (0->0), fcn. (45->8), ass. (0->19)
t48 = qJ(1) + qJ(2);
t45 = sin(t48);
t51 = sin(qJ(1));
t57 = t51 * pkin(1) + pkin(2) * t45;
t56 = qJ(3) + pkin(6) + pkin(5);
t46 = cos(t48);
t53 = cos(qJ(1));
t55 = -pkin(1) * t53 - pkin(2) * t46;
t44 = pkin(8) + t48;
t40 = sin(t44);
t41 = cos(t44);
t54 = g(2) * t40 - g(3) * t41;
t52 = cos(qJ(4));
t50 = sin(qJ(4));
t49 = -qJ(5) - pkin(7);
t43 = pkin(4) * t52 + pkin(3);
t39 = -g(1) * t52 + t50 * t54;
t38 = -g(1) * t50 - t52 * t54;
t1 = [0, -g(2) * t51 + g(3) * t53, -g(2) * t53 - g(3) * t51, 0, -g(2) * t45 + g(3) * t46, -g(2) * t46 - g(3) * t45, -g(1) * t56 - g(2) * t57 - g(3) * t55, 0, 0, 0, 0, 0, t38, t39, t38, t39, g(2) * t41 + g(3) * t40, -g(1) * (pkin(4) * t50 + t56) - g(2) * (t40 * t43 + t41 * t49 + t57) - g(3) * (t40 * t49 - t41 * t43 + t55);];
U_reg = t1;

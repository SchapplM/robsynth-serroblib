% Calculate minimal parameter regressor of potential energy for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:47
% EndTime: 2021-01-15 19:14:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (81->35), mult. (100->46), div. (0->0), fcn. (102->8), ass. (0->24)
t60 = pkin(8) + qJ(3);
t58 = sin(t60);
t76 = g(3) * t58;
t65 = sin(qJ(4));
t66 = sin(qJ(1));
t75 = t66 * t65;
t67 = cos(qJ(4));
t74 = t66 * t67;
t68 = cos(qJ(1));
t73 = t68 * t65;
t72 = t68 * t67;
t71 = pkin(4) * t65 + pkin(6) + qJ(2);
t70 = g(1) * t68 + g(2) * t66;
t57 = t67 * pkin(4) + pkin(3);
t59 = cos(t60);
t62 = cos(pkin(8));
t63 = -qJ(5) - pkin(7);
t69 = t62 * pkin(2) + t57 * t59 - t58 * t63 + pkin(1);
t61 = sin(pkin(8));
t55 = g(1) * t66 - g(2) * t68;
t54 = -g(3) * t59 + t70 * t58;
t53 = -g(1) * (t59 * t72 + t75) - g(2) * (t59 * t74 - t73) - t67 * t76;
t52 = -g(1) * (-t59 * t73 + t74) - g(2) * (-t59 * t75 - t72) + t65 * t76;
t1 = [0, -t70, t55, -g(3) * t61 - t70 * t62, -t55, -g(1) * (t68 * pkin(1) + t66 * qJ(2)) - g(2) * (t66 * pkin(1) - t68 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t70 * t59 - t76, t54, 0, 0, 0, 0, 0, t53, t52, t53, t52, -t54, -g(3) * (t61 * pkin(2) + t58 * t57 + t59 * t63 + pkin(5)) + (-g(1) * t69 + g(2) * t71) * t68 + (-g(1) * t71 - g(2) * t69) * t66;];
U_reg = t1;

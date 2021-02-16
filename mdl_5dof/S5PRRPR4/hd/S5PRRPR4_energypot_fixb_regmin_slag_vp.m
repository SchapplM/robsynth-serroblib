% Calculate minimal parameter regressor of potential energy for
% S5PRRPR4
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
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:40
% EndTime: 2021-01-15 15:52:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (79->42), mult. (100->66), div. (0->0), fcn. (108->10), ass. (0->24)
t62 = sin(qJ(2));
t72 = g(3) * t62;
t58 = sin(pkin(8));
t64 = cos(qJ(2));
t71 = t58 * t64;
t59 = cos(pkin(8));
t70 = t59 * t64;
t61 = sin(qJ(3));
t69 = t61 * t64;
t63 = cos(qJ(3));
t68 = t63 * t64;
t57 = qJ(3) + pkin(9);
t67 = pkin(3) * t61 + pkin(5);
t66 = g(1) * t59 + g(2) * t58;
t53 = t63 * pkin(3) + pkin(2);
t60 = qJ(4) + pkin(6);
t65 = t53 * t64 + t60 * t62 + pkin(1);
t56 = qJ(5) + t57;
t55 = cos(t57);
t54 = sin(t57);
t52 = cos(t56);
t51 = sin(t56);
t50 = -g(3) * t64 + t66 * t62;
t1 = [-g(3) * qJ(1), 0, -t66 * t64 - t72, t50, 0, 0, 0, 0, 0, -g(1) * (t58 * t61 + t59 * t68) - g(2) * (t58 * t68 - t59 * t61) - t63 * t72, -g(1) * (t58 * t63 - t59 * t69) - g(2) * (-t58 * t69 - t59 * t63) + t61 * t72, -g(1) * (t58 * t54 + t55 * t70) - g(2) * (-t59 * t54 + t55 * t71) - t55 * t72, -g(1) * (-t54 * t70 + t58 * t55) - g(2) * (-t54 * t71 - t59 * t55) + t54 * t72, -t50, -g(3) * (t62 * t53 - t64 * t60 + qJ(1)) + (-g(1) * t65 + g(2) * t67) * t59 + (-g(1) * t67 - g(2) * t65) * t58, 0, 0, 0, 0, 0, -g(1) * (t58 * t51 + t52 * t70) - g(2) * (-t59 * t51 + t52 * t71) - t52 * t72, -g(1) * (-t51 * t70 + t58 * t52) - g(2) * (-t51 * t71 - t59 * t52) + t51 * t72;];
U_reg = t1;

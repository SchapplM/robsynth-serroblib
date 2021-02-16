% Calculate minimal parameter regressor of potential energy for
% S5RPRRP3
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:46:04
% EndTime: 2021-01-15 12:46:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (68->23), mult. (56->31), div. (0->0), fcn. (50->8), ass. (0->18)
t66 = qJ(2) + pkin(5);
t58 = qJ(1) + pkin(8);
t53 = sin(t58);
t54 = cos(t58);
t65 = g(2) * t53 - g(3) * t54;
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t64 = -g(2) * t61 + g(3) * t63;
t62 = cos(qJ(3));
t60 = sin(qJ(3));
t59 = qJ(3) + qJ(4);
t57 = qJ(5) + pkin(7) + pkin(6);
t56 = cos(t59);
t55 = sin(t59);
t52 = t62 * pkin(3) + pkin(4) * t56 + pkin(2);
t51 = -g(1) * t56 + t65 * t55;
t50 = -g(1) * t55 - t65 * t56;
t1 = [0, t64, -g(2) * t63 - g(3) * t61, t64 * pkin(1) - g(1) * t66, 0, 0, 0, 0, 0, -g(1) * t60 - t65 * t62, -g(1) * t62 + t65 * t60, 0, 0, 0, 0, 0, t50, t51, t50, t51, g(2) * t54 + g(3) * t53, -g(1) * (t60 * pkin(3) + pkin(4) * t55 + t66) - g(2) * (t61 * pkin(1) + t53 * t52 - t54 * t57) - g(3) * (-t63 * pkin(1) - t54 * t52 - t53 * t57);];
U_reg = t1;

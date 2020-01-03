% Calculate inertial parameters regressor of potential energy for
% S5RRPRP3
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
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:12
% EndTime: 2019-12-31 19:51:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->42), mult. (93->45), div. (0->0), fcn. (80->8), ass. (0->26)
t67 = pkin(6) + pkin(5);
t74 = g(3) * t67;
t62 = sin(pkin(8));
t73 = t62 * pkin(3) + t67;
t63 = cos(pkin(8));
t51 = t63 * pkin(3) + pkin(2);
t61 = qJ(1) + qJ(2);
t55 = sin(t61);
t56 = cos(t61);
t65 = sin(qJ(1));
t58 = t65 * pkin(1);
t64 = -pkin(7) - qJ(3);
t72 = t55 * t51 + t56 * t64 + t58;
t66 = cos(qJ(1));
t59 = t66 * pkin(1);
t71 = t56 * t51 - t55 * t64 + t59;
t70 = g(1) * t56 + g(2) * t55;
t69 = -g(1) * t66 - g(2) * t65;
t60 = pkin(8) + qJ(4);
t53 = sin(t60);
t54 = cos(t60);
t68 = pkin(4) * t54 + qJ(5) * t53;
t47 = g(1) * t55 - g(2) * t56;
t46 = -g(3) * t53 - t70 * t54;
t45 = -g(3) * t54 + t70 * t53;
t1 = [0, 0, 0, 0, 0, 0, t69, g(1) * t65 - g(2) * t66, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t70, t47, -g(3), t69 * pkin(1) - t74, 0, 0, 0, 0, 0, 0, -g(3) * t62 - t70 * t63, -g(3) * t63 + t70 * t62, -t47, -g(1) * (t56 * pkin(2) + t55 * qJ(3) + t59) - g(2) * (t55 * pkin(2) - t56 * qJ(3) + t58) - t74, 0, 0, 0, 0, 0, 0, t46, t45, -t47, -g(1) * t71 - g(2) * t72 - g(3) * t73, 0, 0, 0, 0, 0, 0, t46, -t47, -t45, -g(1) * (t68 * t56 + t71) - g(2) * (t68 * t55 + t72) - g(3) * (t53 * pkin(4) - t54 * qJ(5) + t73);];
U_reg = t1;

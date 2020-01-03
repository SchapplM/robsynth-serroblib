% Calculate inertial parameters regressor of potential energy for
% S5RRPRP4
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
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:55
% EndTime: 2019-12-31 19:52:55
% DurationCPUTime: 0.05s
% Computational Cost: add. (96->38), mult. (83->37), div. (0->0), fcn. (70->6), ass. (0->23)
t55 = pkin(6) + pkin(5);
t50 = qJ(1) + qJ(2);
t47 = cos(t50);
t64 = g(2) * t47;
t63 = g(3) * t55;
t46 = sin(t50);
t52 = sin(qJ(1));
t62 = t52 * pkin(1) + t46 * pkin(2);
t61 = pkin(3) + t55;
t54 = cos(qJ(1));
t60 = t54 * pkin(1) + t47 * pkin(2) + t46 * qJ(3);
t59 = t47 * pkin(7) + t60;
t58 = -t47 * qJ(3) + t62;
t39 = g(1) * t46 - t64;
t57 = -g(1) * t54 - g(2) * t52;
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t56 = pkin(4) * t51 - qJ(5) * t53;
t42 = t46 * pkin(7);
t40 = g(1) * t47 + g(2) * t46;
t38 = -g(3) * t51 + t39 * t53;
t37 = -g(3) * t53 - t39 * t51;
t1 = [0, 0, 0, 0, 0, 0, t57, g(1) * t52 - g(2) * t54, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t40, t39, -g(3), t57 * pkin(1) - t63, 0, 0, 0, 0, 0, 0, -g(3), t40, -t39, -g(1) * t60 - g(2) * t58 - t63, 0, 0, 0, 0, 0, 0, t37, -t38, -t40, -g(1) * t59 - g(2) * (t42 + t58) - g(3) * t61, 0, 0, 0, 0, 0, 0, t37, -t40, t38, -g(1) * (t56 * t46 + t59) - g(2) * (t42 + t62) - g(3) * (t53 * pkin(4) + t51 * qJ(5) + t61) - (-qJ(3) - t56) * t64;];
U_reg = t1;

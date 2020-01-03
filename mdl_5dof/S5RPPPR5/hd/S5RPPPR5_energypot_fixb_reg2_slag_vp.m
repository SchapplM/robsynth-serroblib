% Calculate inertial parameters regressor of potential energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:29
% EndTime: 2019-12-31 17:46:29
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->41), mult. (130->47), div. (0->0), fcn. (141->8), ass. (0->25)
t66 = g(3) * pkin(5);
t53 = -qJ(3) + pkin(5);
t65 = g(3) * t53;
t64 = sin(qJ(1));
t55 = cos(qJ(1));
t63 = t55 * pkin(1) + t64 * qJ(2);
t62 = cos(pkin(7));
t61 = sin(pkin(7));
t60 = t55 * pkin(2) + t63;
t59 = t64 * pkin(1) - t55 * qJ(2);
t58 = t64 * pkin(2) + t59;
t35 = -t55 * t62 - t64 * t61;
t36 = t55 * t61 - t64 * t62;
t57 = g(1) * t36 - g(2) * t35;
t56 = g(1) * t35 + g(2) * t36;
t54 = -pkin(6) - qJ(4);
t52 = cos(pkin(8));
t51 = sin(pkin(8));
t50 = pkin(8) + qJ(5);
t44 = cos(t50);
t43 = sin(t50);
t42 = t52 * pkin(4) + pkin(3);
t38 = -g(1) * t55 - g(2) * t64;
t37 = g(1) * t64 - g(2) * t55;
t1 = [0, 0, 0, 0, 0, 0, t38, t37, -g(3), -t66, 0, 0, 0, 0, 0, 0, t38, -g(3), -t37, -g(1) * t63 - g(2) * t59 - t66, 0, 0, 0, 0, 0, 0, t56, t57, g(3), -g(1) * t60 - g(2) * t58 - t65, 0, 0, 0, 0, 0, 0, g(3) * t51 + t56 * t52, g(3) * t52 - t56 * t51, -t57, -g(1) * (-t35 * pkin(3) + t36 * qJ(4) + t60) - g(2) * (-t36 * pkin(3) - t35 * qJ(4) + t58) - t65, 0, 0, 0, 0, 0, 0, g(3) * t43 + t56 * t44, g(3) * t44 - t56 * t43, -t57, -g(1) * (-t35 * t42 - t36 * t54 + t60) - g(2) * (t35 * t54 - t36 * t42 + t58) - g(3) * (-t51 * pkin(4) + t53);];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:38
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->36), mult. (97->42), div. (0->0), fcn. (98->8), ass. (0->24)
t45 = qJ(2) + pkin(5);
t62 = g(3) * (-pkin(6) + t45);
t61 = g(3) * t45;
t60 = cos(qJ(4));
t59 = sin(qJ(4));
t58 = qJ(1) + pkin(8);
t41 = cos(t58);
t49 = cos(qJ(1));
t55 = sin(t58);
t57 = t49 * pkin(1) + t41 * pkin(2) + t55 * qJ(3);
t56 = t41 * pkin(3) + t57;
t47 = sin(qJ(1));
t54 = t47 * pkin(1) + t55 * pkin(2) - t41 * qJ(3);
t29 = -t41 * t60 - t55 * t59;
t30 = t41 * t59 - t55 * t60;
t53 = g(1) * t30 - g(2) * t29;
t52 = g(1) * t29 + g(2) * t30;
t51 = -g(1) * t49 - g(2) * t47;
t50 = t55 * pkin(3) + t54;
t48 = cos(qJ(5));
t46 = sin(qJ(5));
t32 = -g(1) * t41 - g(2) * t55;
t31 = g(1) * t55 - g(2) * t41;
t1 = [0, 0, 0, 0, 0, 0, t51, g(1) * t47 - g(2) * t49, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, t32, t31, -g(3), t51 * pkin(1) - t61, 0, 0, 0, 0, 0, 0, t32, -g(3), -t31, -g(1) * t57 - g(2) * t54 - t61, 0, 0, 0, 0, 0, 0, t52, t53, g(3), -g(1) * t56 - g(2) * t50 - t62, 0, 0, 0, 0, 0, 0, g(3) * t46 + t52 * t48, g(3) * t48 - t52 * t46, -t53, -g(1) * (-t29 * pkin(4) + t30 * pkin(7) + t56) - g(2) * (-t30 * pkin(4) - t29 * pkin(7) + t50) - t62;];
U_reg = t1;

% Calculate inertial parameters regressor of potential energy for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:48
% EndTime: 2019-12-31 19:27:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->36), mult. (97->42), div. (0->0), fcn. (98->8), ass. (0->24)
t52 = pkin(6) + pkin(5);
t65 = g(3) * (-qJ(4) + t52);
t64 = g(3) * t52;
t63 = cos(pkin(8));
t62 = sin(pkin(8));
t61 = qJ(1) + qJ(2);
t44 = cos(t61);
t51 = cos(qJ(1));
t59 = sin(t61);
t60 = t51 * pkin(1) + t44 * pkin(2) + t59 * qJ(3);
t58 = t44 * pkin(3) + t60;
t49 = sin(qJ(1));
t57 = t49 * pkin(1) + t59 * pkin(2) - t44 * qJ(3);
t32 = -t44 * t63 - t59 * t62;
t33 = t44 * t62 - t59 * t63;
t56 = g(1) * t33 - g(2) * t32;
t55 = g(1) * t32 + g(2) * t33;
t54 = -g(1) * t51 - g(2) * t49;
t53 = t59 * pkin(3) + t57;
t50 = cos(qJ(5));
t48 = sin(qJ(5));
t35 = -g(1) * t44 - g(2) * t59;
t34 = g(1) * t59 - g(2) * t44;
t1 = [0, 0, 0, 0, 0, 0, t54, g(1) * t49 - g(2) * t51, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, t35, t34, -g(3), t54 * pkin(1) - t64, 0, 0, 0, 0, 0, 0, t35, -g(3), -t34, -g(1) * t60 - g(2) * t57 - t64, 0, 0, 0, 0, 0, 0, t55, t56, g(3), -g(1) * t58 - g(2) * t53 - t65, 0, 0, 0, 0, 0, 0, g(3) * t48 + t55 * t50, g(3) * t50 - t55 * t48, -t56, -g(1) * (-t32 * pkin(4) + t33 * pkin(7) + t58) - g(2) * (-t33 * pkin(4) - t32 * pkin(7) + t53) - t65;];
U_reg = t1;

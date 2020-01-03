% Calculate inertial parameters regressor of potential energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:11
% EndTime: 2019-12-31 17:41:11
% DurationCPUTime: 0.06s
% Computational Cost: add. (109->40), mult. (104->43), div. (0->0), fcn. (91->6), ass. (0->24)
t49 = pkin(7) + qJ(2);
t43 = sin(t49);
t53 = sin(qJ(3));
t61 = qJ(4) * t53;
t54 = cos(qJ(3));
t64 = t43 * t54;
t66 = pkin(3) * t64 + t43 * t61;
t52 = pkin(5) + qJ(1);
t65 = g(3) * t52;
t44 = cos(t49);
t63 = t44 * t54;
t50 = sin(pkin(7));
t62 = t50 * pkin(1) + t43 * pkin(2);
t51 = cos(pkin(7));
t60 = t51 * pkin(1) + t44 * pkin(2) + t43 * pkin(6);
t59 = -t44 * pkin(6) + t62;
t58 = pkin(3) * t63 + t44 * t61 + t60;
t57 = g(1) * t44 + g(2) * t43;
t56 = -g(1) * t51 - g(2) * t50;
t55 = t53 * pkin(3) - t54 * qJ(4) + t52;
t34 = g(1) * t43 - g(2) * t44;
t33 = -g(3) * t53 - t57 * t54;
t32 = -g(3) * t54 + t57 * t53;
t1 = [0, 0, 0, 0, 0, 0, t56, g(1) * t50 - g(2) * t51, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t57, t34, -g(3), t56 * pkin(1) - t65, 0, 0, 0, 0, 0, 0, t33, t32, -t34, -g(1) * t60 - g(2) * t59 - t65, 0, 0, 0, 0, 0, 0, t33, -t34, -t32, -g(1) * t58 - g(2) * (t59 + t66) - g(3) * t55, 0, 0, 0, 0, 0, 0, t33, -t32, t34, -g(1) * (pkin(4) * t63 - t43 * qJ(5) + t58) - g(2) * (pkin(4) * t64 + (-pkin(6) + qJ(5)) * t44 + t62 + t66) - g(3) * (t53 * pkin(4) + t55);];
U_reg = t1;

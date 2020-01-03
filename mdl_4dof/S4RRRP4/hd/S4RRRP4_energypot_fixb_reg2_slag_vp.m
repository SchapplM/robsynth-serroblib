% Calculate inertial parameters regressor of potential energy for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->31), mult. (71->37), div. (0->0), fcn. (61->6), ass. (0->18)
t58 = g(3) * pkin(4);
t55 = -pkin(6) - pkin(5);
t51 = sin(qJ(2));
t57 = t51 * pkin(2) + pkin(4);
t53 = cos(qJ(2));
t43 = t53 * pkin(2) + pkin(1);
t52 = sin(qJ(1));
t54 = cos(qJ(1));
t56 = g(1) * t54 + g(2) * t52;
t50 = qJ(2) + qJ(3);
t49 = -qJ(4) + t55;
t45 = cos(t50);
t44 = sin(t50);
t42 = g(1) * t52 - g(2) * t54;
t41 = pkin(3) * t45 + t43;
t40 = -g(3) * t44 - t56 * t45;
t39 = -g(3) * t45 + t56 * t44;
t1 = [0, 0, 0, 0, 0, 0, -t56, t42, -g(3), -t58, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t56 * t53, -g(3) * t53 + t56 * t51, -t42, -g(1) * (t54 * pkin(1) + t52 * pkin(5)) - g(2) * (t52 * pkin(1) - t54 * pkin(5)) - t58, 0, 0, 0, 0, 0, 0, t40, t39, -t42, -g(1) * (t54 * t43 - t52 * t55) - g(2) * (t52 * t43 + t54 * t55) - g(3) * t57, 0, 0, 0, 0, 0, 0, t40, t39, -t42, -g(1) * (t54 * t41 - t52 * t49) - g(2) * (t52 * t41 + t54 * t49) - g(3) * (pkin(3) * t44 + t57);];
U_reg = t1;

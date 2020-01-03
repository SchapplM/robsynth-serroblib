% Calculate inertial parameters regressor of potential energy for
% S4RRRP5
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:05
% DurationCPUTime: 0.05s
% Computational Cost: add. (65->32), mult. (78->37), div. (0->0), fcn. (68->6), ass. (0->19)
t60 = g(3) * pkin(4);
t50 = sin(qJ(2));
t59 = t50 * pkin(2) + pkin(4);
t52 = cos(qJ(2));
t43 = t52 * pkin(2) + pkin(1);
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t54 = -pkin(6) - pkin(5);
t58 = t51 * t43 + t53 * t54;
t57 = t53 * t43 - t51 * t54;
t56 = g(1) * t53 + g(2) * t51;
t49 = qJ(2) + qJ(3);
t45 = sin(t49);
t46 = cos(t49);
t55 = pkin(3) * t46 + qJ(4) * t45;
t40 = g(1) * t51 - g(2) * t53;
t39 = -g(3) * t45 - t56 * t46;
t38 = -g(3) * t46 + t56 * t45;
t1 = [0, 0, 0, 0, 0, 0, -t56, t40, -g(3), -t60, 0, 0, 0, 0, 0, 0, -g(3) * t50 - t56 * t52, -g(3) * t52 + t56 * t50, -t40, -g(1) * (t53 * pkin(1) + t51 * pkin(5)) - g(2) * (t51 * pkin(1) - t53 * pkin(5)) - t60, 0, 0, 0, 0, 0, 0, t39, t38, -t40, -g(1) * t57 - g(2) * t58 - g(3) * t59, 0, 0, 0, 0, 0, 0, t39, -t40, -t38, -g(1) * (t55 * t53 + t57) - g(2) * (t55 * t51 + t58) - g(3) * (t45 * pkin(3) - t46 * qJ(4) + t59);];
U_reg = t1;

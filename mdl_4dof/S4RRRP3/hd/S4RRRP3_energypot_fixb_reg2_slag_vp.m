% Calculate inertial parameters regressor of potential energy for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:14
% EndTime: 2019-12-31 17:14:14
% DurationCPUTime: 0.07s
% Computational Cost: add. (68->29), mult. (66->32), div. (0->0), fcn. (56->6), ass. (0->18)
t48 = pkin(5) + pkin(4);
t54 = g(3) * t48;
t43 = qJ(1) + qJ(2);
t39 = sin(t43);
t40 = cos(t43);
t47 = cos(qJ(1));
t53 = t47 * pkin(1) + t40 * pkin(2) + t39 * pkin(6);
t45 = sin(qJ(1));
t52 = t45 * pkin(1) + t39 * pkin(2) - t40 * pkin(6);
t51 = g(1) * t40 + g(2) * t39;
t50 = -g(1) * t47 - g(2) * t45;
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t49 = pkin(3) * t46 + qJ(4) * t44;
t34 = g(1) * t39 - g(2) * t40;
t33 = -g(3) * t44 - t51 * t46;
t32 = -g(3) * t46 + t51 * t44;
t1 = [0, 0, 0, 0, 0, 0, t50, g(1) * t45 - g(2) * t47, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t51, t34, -g(3), t50 * pkin(1) - t54, 0, 0, 0, 0, 0, 0, t33, t32, -t34, -g(1) * t53 - g(2) * t52 - t54, 0, 0, 0, 0, 0, 0, t33, -t34, -t32, -g(1) * (t49 * t40 + t53) - g(2) * (t49 * t39 + t52) - g(3) * (t44 * pkin(3) - t46 * qJ(4) + t48);];
U_reg = t1;

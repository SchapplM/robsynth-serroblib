% Calculate inertial parameters regressor of potential energy for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:17
% EndTime: 2019-12-31 16:59:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (48->32), mult. (89->35), div. (0->0), fcn. (79->4), ass. (0->19)
t45 = sin(qJ(1));
t44 = sin(qJ(2));
t52 = qJ(3) * t44;
t46 = cos(qJ(2));
t55 = t45 * t46;
t57 = pkin(2) * t55 + t45 * t52;
t56 = g(3) * pkin(4);
t47 = cos(qJ(1));
t54 = t46 * t47;
t53 = t47 * pkin(1) + t45 * pkin(5);
t40 = t45 * pkin(1);
t51 = -t47 * pkin(5) + t40;
t50 = pkin(2) * t54 + t47 * t52 + t53;
t49 = t44 * pkin(2) - t46 * qJ(3) + pkin(4);
t48 = g(1) * t47 + g(2) * t45;
t33 = g(1) * t45 - g(2) * t47;
t32 = -g(3) * t44 - t48 * t46;
t31 = -g(3) * t46 + t44 * t48;
t1 = [0, 0, 0, 0, 0, 0, -t48, t33, -g(3), -t56, 0, 0, 0, 0, 0, 0, t32, t31, -t33, -g(1) * t53 - g(2) * t51 - t56, 0, 0, 0, 0, 0, 0, t32, -t33, -t31, -g(1) * t50 - g(2) * (t51 + t57) - g(3) * t49, 0, 0, 0, 0, 0, 0, t32, -t31, t33, -g(1) * (pkin(3) * t54 - t45 * qJ(4) + t50) - g(2) * (pkin(3) * t55 + t40 + (-pkin(5) + qJ(4)) * t47 + t57) - g(3) * (t44 * pkin(3) + t49);];
U_reg = t1;

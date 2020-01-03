% Calculate inertial parameters regressor of potential energy for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:36
% EndTime: 2019-12-31 17:00:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (48->34), mult. (89->33), div. (0->0), fcn. (79->4), ass. (0->18)
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t56 = pkin(2) * t44 + qJ(3) * t42;
t43 = sin(qJ(1));
t55 = t56 * t43;
t54 = g(3) * pkin(4);
t45 = cos(qJ(1));
t52 = t45 * pkin(1) + t43 * pkin(5);
t50 = qJ(4) * t44;
t38 = t43 * pkin(1);
t49 = -t45 * pkin(5) + t38;
t48 = t56 * t45 + t52;
t47 = t42 * pkin(2) - t44 * qJ(3) + pkin(4);
t46 = g(1) * t45 + g(2) * t43;
t31 = g(1) * t43 - g(2) * t45;
t30 = g(3) * t42 + t46 * t44;
t29 = -g(3) * t44 + t46 * t42;
t1 = [0, 0, 0, 0, 0, 0, -t46, t31, -g(3), -t54, 0, 0, 0, 0, 0, 0, -t30, t29, -t31, -g(1) * t52 - g(2) * t49 - t54, 0, 0, 0, 0, 0, 0, -t31, t30, -t29, -g(1) * t48 - g(2) * (t49 + t55) - g(3) * t47, 0, 0, 0, 0, 0, 0, -t31, -t29, -t30, -g(1) * (t43 * pkin(3) + t45 * t50 + t48) - g(2) * (t43 * t50 + t38 + (-pkin(3) - pkin(5)) * t45 + t55) - g(3) * (t42 * qJ(4) + t47);];
U_reg = t1;

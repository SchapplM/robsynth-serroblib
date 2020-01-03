% Calculate inertial parameters regressor of potential energy for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:39
% EndTime: 2019-12-31 16:51:39
% DurationCPUTime: 0.05s
% Computational Cost: add. (50->28), mult. (82->34), div. (0->0), fcn. (86->6), ass. (0->19)
t52 = g(3) * pkin(4);
t51 = g(3) * (-pkin(5) + pkin(4));
t50 = cos(qJ(3));
t49 = sin(qJ(1));
t48 = sin(qJ(3));
t40 = cos(qJ(1));
t47 = t40 * pkin(1) + t49 * qJ(2);
t46 = t40 * pkin(2) + t47;
t45 = t49 * pkin(1) - t40 * qJ(2);
t44 = t49 * pkin(2) + t45;
t26 = -t40 * t50 - t49 * t48;
t27 = t40 * t48 - t49 * t50;
t43 = g(1) * t27 - g(2) * t26;
t42 = g(1) * t26 + g(2) * t27;
t39 = cos(qJ(4));
t38 = sin(qJ(4));
t29 = -g(1) * t40 - g(2) * t49;
t28 = g(1) * t49 - g(2) * t40;
t1 = [0, 0, 0, 0, 0, 0, t29, t28, -g(3), -t52, 0, 0, 0, 0, 0, 0, t29, -g(3), -t28, -g(1) * t47 - g(2) * t45 - t52, 0, 0, 0, 0, 0, 0, t42, t43, g(3), -g(1) * t46 - g(2) * t44 - t51, 0, 0, 0, 0, 0, 0, g(3) * t38 + t42 * t39, g(3) * t39 - t42 * t38, -t43, -g(1) * (-t26 * pkin(3) + t27 * pkin(6) + t46) - g(2) * (-t27 * pkin(3) - t26 * pkin(6) + t44) - t51;];
U_reg = t1;

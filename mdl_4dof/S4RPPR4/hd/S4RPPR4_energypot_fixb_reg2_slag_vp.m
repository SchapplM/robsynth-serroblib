% Calculate inertial parameters regressor of potential energy for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:53
% EndTime: 2019-12-31 16:38:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->27), mult. (50->28), div. (0->0), fcn. (40->6), ass. (0->15)
t30 = qJ(2) + pkin(4);
t38 = g(3) * t30;
t29 = qJ(1) + pkin(6);
t25 = sin(t29);
t26 = cos(t29);
t34 = cos(qJ(1));
t37 = t34 * pkin(1) + t26 * pkin(2) + t25 * qJ(3);
t32 = sin(qJ(1));
t36 = t32 * pkin(1) + t25 * pkin(2) - t26 * qJ(3);
t20 = g(1) * t25 - g(2) * t26;
t35 = -g(1) * t34 - g(2) * t32;
t33 = cos(qJ(4));
t31 = sin(qJ(4));
t21 = g(1) * t26 + g(2) * t25;
t1 = [0, 0, 0, 0, 0, 0, t35, g(1) * t32 - g(2) * t34, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t21, t20, -g(3), t35 * pkin(1) - t38, 0, 0, 0, 0, 0, 0, -g(3), t21, -t20, -g(1) * t37 - g(2) * t36 - t38, 0, 0, 0, 0, 0, 0, -g(3) * t33 - t20 * t31, g(3) * t31 - t20 * t33, -t21, -g(1) * (t26 * pkin(5) + t37) - g(2) * (t25 * pkin(5) + t36) - g(3) * (pkin(3) + t30);];
U_reg = t1;

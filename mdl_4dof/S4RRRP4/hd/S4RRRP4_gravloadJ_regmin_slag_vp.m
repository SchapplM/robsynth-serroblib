% Calculate minimal parameter regressor of gravitation load for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:22
% EndTime: 2021-01-15 14:30:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (88->20), mult. (111->29), div. (0->0), fcn. (103->6), ass. (0->17)
t14 = qJ(2) + qJ(3);
t11 = cos(t14);
t17 = cos(qJ(2));
t19 = t17 * pkin(2) + pkin(3) * t11;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t8 = g(1) * t18 + g(2) * t16;
t7 = g(1) * t16 - g(2) * t18;
t10 = sin(t14);
t1 = -g(3) * t11 + t8 * t10;
t15 = sin(qJ(2));
t13 = -qJ(4) - pkin(6) - pkin(5);
t5 = pkin(1) + t19;
t4 = t7 * t11;
t3 = t7 * t10;
t2 = g(3) * t10 + t8 * t11;
t6 = [0, t7, t8, 0, 0, 0, 0, 0, t7 * t17, -t7 * t15, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t8, -g(1) * (-t18 * t13 - t16 * t5) - g(2) * (-t16 * t13 + t18 * t5); 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t17 + t8 * t15, g(3) * t15 + t8 * t17, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(3) * t19 - t8 * (-t15 * pkin(2) - pkin(3) * t10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t6;

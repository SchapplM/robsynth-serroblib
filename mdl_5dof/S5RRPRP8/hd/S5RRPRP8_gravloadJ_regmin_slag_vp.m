% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:19
% EndTime: 2021-01-15 20:52:20
% DurationCPUTime: 0.13s
% Computational Cost: add. (121->38), mult. (245->47), div. (0->0), fcn. (264->6), ass. (0->26)
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t36 = sin(qJ(4));
t37 = cos(qJ(4));
t14 = -t30 * t36 - t32 * t37;
t15 = t30 * t37 - t32 * t36;
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t20 = g(1) * t33 + g(2) * t31;
t3 = g(3) * t14 + t20 * t15;
t4 = -g(3) * t15 + t20 * t14;
t35 = t32 * pkin(2) + t30 * qJ(3);
t21 = t37 * pkin(4) + pkin(2) + pkin(3);
t22 = t36 * pkin(4) + qJ(3);
t34 = t21 * t32 + t22 * t30;
t19 = g(1) * t31 - g(2) * t33;
t29 = qJ(5) - pkin(6) + pkin(7);
t16 = pkin(1) + t35;
t13 = t19 * t32;
t12 = t19 * t30;
t10 = g(3) * t30 + t20 * t32;
t9 = -g(3) * t32 + t20 * t30;
t7 = pkin(1) + t34;
t6 = t19 * t14;
t5 = t19 * t15;
t1 = [0, t19, t20, 0, 0, 0, 0, 0, t13, -t12, t13, -t20, t12, -g(1) * (t33 * pkin(6) - t16 * t31) - g(2) * (t31 * pkin(6) + t16 * t33), 0, 0, 0, 0, 0, -t6, t5, -t6, t5, t20, -g(1) * (-t29 * t33 - t7 * t31) - g(2) * (-t29 * t31 + t7 * t33); 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, t9, 0, -t10, -g(3) * t35 - t20 * (-t30 * pkin(2) + t32 * qJ(3)), 0, 0, 0, 0, 0, t3, t4, t3, t4, 0, -g(3) * t34 - t20 * (-t21 * t30 + t22 * t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -t3, -t4, 0, -t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19;];
taug_reg = t1;

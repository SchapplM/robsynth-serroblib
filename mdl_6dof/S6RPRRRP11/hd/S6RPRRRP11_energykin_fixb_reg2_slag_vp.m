% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:11
% EndTime: 2019-03-09 06:41:11
% DurationCPUTime: 0.32s
% Computational Cost: add. (1680->78), mult. (5438->177), div. (0->0), fcn. (4485->12), ass. (0->65)
t62 = sin(pkin(12));
t65 = cos(pkin(12));
t64 = sin(pkin(6));
t80 = qJD(1) * t64;
t74 = qJ(2) * t80;
t67 = cos(pkin(6));
t79 = qJD(1) * t67;
t77 = pkin(1) * t79;
t53 = t62 * t77 + t65 * t74;
t63 = sin(pkin(7));
t66 = cos(pkin(7));
t82 = t64 * t65;
t40 = (t63 * t67 + t66 * t82) * qJD(1) * pkin(9) + t53;
t58 = t65 * t77;
t84 = t62 * t64;
t42 = t58 + (pkin(2) * t67 + (-pkin(9) * t66 - qJ(2)) * t84) * qJD(1);
t48 = qJD(2) + (-pkin(9) * t62 * t63 - pkin(2) * t65 - pkin(1)) * t80;
t70 = sin(qJ(3));
t71 = cos(qJ(3));
t23 = -t70 * t40 + (t42 * t66 + t48 * t63) * t71;
t72 = qJD(1) ^ 2;
t88 = t72 / 0.2e1;
t76 = t65 * t80;
t50 = t63 * t76 - t66 * t79 - qJD(3);
t21 = t50 * pkin(3) - t23;
t81 = t66 * t70;
t83 = t63 * t70;
t45 = (t67 * t83 + (t62 * t71 + t65 * t81) * t64) * qJD(1);
t69 = sin(qJ(4));
t87 = cos(qJ(4));
t33 = t69 * t45 + t87 * t50;
t35 = t87 * t45 - t69 * t50;
t13 = t33 * pkin(4) - t35 * pkin(11) + t21;
t68 = sin(qJ(5));
t30 = -t63 * t42 + t66 * t48;
t43 = t70 * t62 * t80 + (-t63 * t79 - t66 * t76) * t71;
t18 = t43 * pkin(3) - t45 * pkin(10) + t30;
t24 = t71 * t40 + t42 * t81 + t48 * t83;
t22 = -t50 * pkin(10) + t24;
t10 = t69 * t18 + t87 * t22;
t41 = qJD(4) + t43;
t8 = t41 * pkin(11) + t10;
t86 = cos(qJ(5));
t4 = t68 * t13 + t86 * t8;
t85 = t64 ^ 2 * t72;
t78 = t64 * t67 * t72;
t75 = t85 / 0.2e1;
t3 = t86 * t13 - t68 * t8;
t9 = t87 * t18 - t69 * t22;
t7 = -t41 * pkin(4) - t9;
t59 = -pkin(1) * t80 + qJD(2);
t52 = -t62 * t74 + t58;
t32 = qJD(5) + t33;
t31 = t32 ^ 2 / 0.2e1;
t29 = t86 * t35 + t68 * t41;
t27 = t68 * t35 - t86 * t41;
t26 = t29 ^ 2 / 0.2e1;
t25 = t27 ^ 2 / 0.2e1;
t16 = t29 * t32;
t15 = t27 * t32;
t14 = t29 * t27;
t5 = t27 * pkin(5) + qJD(6) + t7;
t2 = -t27 * qJ(6) + t4;
t1 = t32 * pkin(5) - t29 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t88, 0, 0, 0, 0, t62 ^ 2 * t75, t62 * t65 * t85, t62 * t78, t65 ^ 2 * t75, t65 * t78, t67 ^ 2 * t88 (t52 * t67 - t59 * t82) * qJD(1) (-t53 * t67 + t59 * t84) * qJD(1) (-t52 * t62 + t53 * t65) * t80, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t23 * t50 + t30 * t43, t24 * t50 + t30 * t45, -t23 * t45 - t24 * t43, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t41, t33 ^ 2 / 0.2e1, -t33 * t41, t41 ^ 2 / 0.2e1, t21 * t33 + t9 * t41, -t10 * t41 + t21 * t35, -t10 * t33 - t9 * t35, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t26, -t14, t16, t25, -t15, t31, t7 * t27 + t3 * t32, t7 * t29 - t4 * t32, -t4 * t27 - t3 * t29, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t26, -t14, t16, t25, -t15, t31, t1 * t32 + t5 * t27, -t2 * t32 + t5 * t29, -t1 * t29 - t2 * t27, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;

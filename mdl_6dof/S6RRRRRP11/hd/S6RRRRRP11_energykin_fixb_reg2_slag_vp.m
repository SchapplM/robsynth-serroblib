% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP11_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energykin_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:35
% EndTime: 2019-03-10 02:53:36
% DurationCPUTime: 0.34s
% Computational Cost: add. (2033->76), mult. (5439->166), div. (0->0), fcn. (4483->12), ass. (0->62)
t80 = cos(pkin(6)) * qJD(1);
t60 = qJD(2) + t80;
t62 = sin(pkin(7));
t64 = cos(pkin(7));
t71 = cos(qJ(2));
t63 = sin(pkin(6));
t81 = qJD(1) * t63;
t76 = t71 * t81;
t88 = t60 * t62 + t64 * t76;
t69 = sin(qJ(2));
t78 = pkin(1) * t80;
t53 = pkin(9) * t76 + t69 * t78;
t40 = t88 * pkin(10) + t53;
t59 = t71 * t78;
t77 = t69 * t81;
t42 = t60 * pkin(2) + t59 + (-pkin(10) * t64 - pkin(9)) * t77;
t49 = (-pkin(10) * t62 * t69 - pkin(2) * t71 - pkin(1)) * t81;
t68 = sin(qJ(3));
t70 = cos(qJ(3));
t23 = -t68 * t40 + (t42 * t64 + t49 * t62) * t70;
t50 = -t64 * t60 + t62 * t76 - qJD(3);
t21 = t50 * pkin(3) - t23;
t82 = t64 * t68;
t83 = t62 * t68;
t45 = t60 * t83 + (t69 * t70 + t71 * t82) * t81;
t67 = sin(qJ(4));
t87 = cos(qJ(4));
t33 = t67 * t45 + t87 * t50;
t35 = t45 * t87 - t67 * t50;
t13 = t33 * pkin(4) - t35 * pkin(12) + t21;
t66 = sin(qJ(5));
t30 = -t62 * t42 + t64 * t49;
t43 = t68 * t77 - t88 * t70;
t18 = t43 * pkin(3) - t45 * pkin(11) + t30;
t24 = t70 * t40 + t42 * t82 + t49 * t83;
t22 = -t50 * pkin(11) + t24;
t10 = t67 * t18 + t87 * t22;
t41 = qJD(4) + t43;
t8 = t41 * pkin(12) + t10;
t86 = cos(qJ(5));
t4 = t66 * t13 + t86 * t8;
t72 = qJD(1) ^ 2;
t84 = t63 ^ 2 * t72;
t79 = t71 * t84;
t75 = t84 / 0.2e1;
t3 = t86 * t13 - t66 * t8;
t9 = t18 * t87 - t67 * t22;
t7 = -t41 * pkin(4) - t9;
t52 = -pkin(9) * t77 + t59;
t32 = qJD(5) + t33;
t31 = t32 ^ 2 / 0.2e1;
t29 = t35 * t86 + t66 * t41;
t27 = t66 * t35 - t41 * t86;
t26 = t29 ^ 2 / 0.2e1;
t25 = t27 ^ 2 / 0.2e1;
t16 = t29 * t32;
t15 = t27 * t32;
t14 = t29 * t27;
t5 = t27 * pkin(5) + qJD(6) + t7;
t2 = -t27 * qJ(6) + t4;
t1 = t32 * pkin(5) - t29 * qJ(6) + t3;
t6 = [0, 0, 0, 0, 0, t72 / 0.2e1, 0, 0, 0, 0, t69 ^ 2 * t75, t69 * t79, t60 * t77, t71 ^ 2 * t75, t60 * t76, t60 ^ 2 / 0.2e1, pkin(1) * t79 + t52 * t60, -pkin(1) * t69 * t84 - t53 * t60 (-t52 * t69 + t53 * t71) * t81, t53 ^ 2 / 0.2e1 + t52 ^ 2 / 0.2e1 + pkin(1) ^ 2 * t75, t45 ^ 2 / 0.2e1, -t45 * t43, -t45 * t50, t43 ^ 2 / 0.2e1, t43 * t50, t50 ^ 2 / 0.2e1, -t23 * t50 + t30 * t43, t24 * t50 + t30 * t45, -t23 * t45 - t24 * t43, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t33, t35 * t41, t33 ^ 2 / 0.2e1, -t33 * t41, t41 ^ 2 / 0.2e1, t21 * t33 + t9 * t41, -t10 * t41 + t21 * t35, -t10 * t33 - t9 * t35, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1, t26, -t14, t16, t25, -t15, t31, t7 * t27 + t3 * t32, t7 * t29 - t4 * t32, -t4 * t27 - t3 * t29, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t26, -t14, t16, t25, -t15, t31, t1 * t32 + t5 * t27, -t2 * t32 + t5 * t29, -t1 * t29 - t2 * t27, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t6;

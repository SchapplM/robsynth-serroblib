% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:04
% EndTime: 2019-12-31 16:45:05
% DurationCPUTime: 0.53s
% Computational Cost: add. (564->66), mult. (1156->78), div. (0->0), fcn. (1262->4), ass. (0->51)
t46 = sin(pkin(6));
t47 = cos(pkin(6));
t48 = sin(qJ(3));
t72 = cos(qJ(3));
t35 = t48 * t46 - t72 * t47;
t31 = t35 ^ 2;
t37 = t72 * t46 + t48 * t47;
t73 = t37 ^ 2;
t10 = t73 - t31;
t80 = t10 * qJD(1);
t79 = t10 * qJD(3);
t76 = t31 + t73;
t78 = t76 * qJD(1);
t77 = t76 * qJD(2);
t71 = pkin(5) + qJ(2);
t40 = t71 * t47;
t51 = t71 * t46;
t20 = t48 * t40 + t72 * t51;
t21 = t72 * t40 - t48 * t51;
t49 = t20 * t37 - t21 * t35;
t75 = t49 * qJD(1);
t74 = t49 * qJD(2);
t41 = t46 ^ 2 + t47 ^ 2;
t43 = -t47 * pkin(2) - pkin(1);
t50 = t35 * pkin(3) - t37 * qJ(4);
t12 = t43 + t50;
t14 = pkin(3) * t37 + t35 * qJ(4);
t1 = t12 * t14;
t70 = t1 * qJD(1);
t5 = t12 * t37 + t14 * t35;
t69 = t5 * qJD(1);
t6 = t12 * t35 - t14 * t37;
t68 = t6 * qJD(1);
t65 = t14 * qJD(1);
t60 = t20 * qJD(3);
t15 = t21 * qJD(3);
t59 = t73 * qJD(1);
t58 = t35 * qJD(1);
t25 = t35 * qJD(3);
t57 = t37 * qJD(1);
t27 = t37 * qJD(3);
t56 = t37 * qJD(4);
t39 = t41 * qJ(2);
t55 = t39 * qJD(1);
t54 = t41 * qJD(1);
t53 = qJD(3) * qJ(4);
t19 = t35 * t57;
t18 = t35 * t27;
t52 = t43 * t57;
t29 = t37 * qJD(2);
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * qJD(2), t39 * qJD(2), -t18, -t79, 0, t18, 0, 0, t43 * t27, -t43 * t25, t77, t74, -t18, 0, t79, 0, 0, t18, t5 * qJD(3) - t35 * t56, t77, t6 * qJD(3) + qJD(4) * t73, t1 * qJD(3) - t12 * t56 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t55, 0, 0, 0, 0, 0, 0, 0, 0, t78, t75, 0, 0, 0, 0, 0, 0, 0, t78, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t80, -t25, t19, -t27, 0, -t15 + t52, -t43 * t58 + t60, 0, 0, -t19, -t25, t80, 0, t27, t19, -t15 + t69, t50 * qJD(3) - t35 * qJD(4), -t60 + t68, t70 + (-t21 * pkin(3) - t20 * qJ(4)) * qJD(3) + t21 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t25, t59, -t12 * t57 + t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0, 0, 0, 0, 0, 0, t27, -t25, -t78, -t75, 0, 0, 0, 0, 0, 0, t27, -t78, t25, qJD(3) * t14 - t56 - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t58, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t58, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t80, 0, -t19, 0, 0, -t29 - t52, (qJD(1) * t43 + qJD(2)) * t35, 0, 0, t19, 0, -t80, 0, 0, -t19, -t29 - t69, 0, -t35 * qJD(2) - t68, -qJD(2) * t14 - t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t58, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, -t58, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t59, (qJD(1) * t12 + qJD(2)) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;

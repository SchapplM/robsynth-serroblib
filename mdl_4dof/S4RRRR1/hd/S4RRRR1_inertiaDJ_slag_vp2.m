% Calculate time derivative of joint inertia matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:09
% DurationCPUTime: 0.30s
% Computational Cost: add. (309->71), mult. (859->118), div. (0->0), fcn. (487->6), ass. (0->53)
t30 = sin(qJ(4));
t28 = t30 ^ 2;
t33 = cos(qJ(4));
t29 = t33 ^ 2;
t63 = t28 + t29;
t62 = Ifges(5,1) - Ifges(5,2);
t34 = cos(qJ(3));
t61 = t63 * t34;
t18 = -t33 * mrSges(5,1) + t30 * mrSges(5,2);
t31 = sin(qJ(3));
t47 = qJD(3) * t31;
t15 = pkin(2) * t18 * t47;
t46 = qJD(3) * t34;
t43 = pkin(2) * t46;
t41 = mrSges(5,3) * t43;
t20 = t28 * t41;
t21 = t29 * t41;
t39 = mrSges(5,1) * t30 + mrSges(5,2) * t33;
t17 = t39 * qJD(4);
t25 = -t34 * pkin(2) - pkin(3);
t8 = t25 * t17;
t60 = t15 + t20 + t21 + t8;
t59 = pkin(3) * t17;
t35 = cos(qJ(2));
t26 = t35 * pkin(1) + pkin(2);
t32 = sin(qJ(2));
t50 = t31 * t32;
t5 = t26 * t46 + (-t32 * t47 + (t34 * t35 - t50) * qJD(2)) * pkin(1);
t58 = t28 * t5;
t57 = t29 * t5;
t56 = t5 * mrSges(4,2);
t53 = Ifges(5,6) * t30;
t49 = t32 * t34;
t12 = pkin(1) * t49 + t31 * t26;
t48 = pkin(2) * qJD(3);
t44 = qJD(4) * t33;
t42 = t63 * t5;
t40 = -t31 * mrSges(4,1) - t34 * mrSges(4,2);
t11 = -pkin(1) * t50 + t34 * t26;
t38 = (-0.2e1 * Ifges(5,4) * t30 + t62 * t33) * qJD(4) * t30 + (0.2e1 * Ifges(5,4) * t33 + t62 * t30) * t44;
t37 = (-mrSges(3,1) * t32 - mrSges(3,2) * t35) * qJD(2) * pkin(1);
t6 = t26 * t47 + (t32 * t46 + (t31 * t35 + t49) * qJD(2)) * pkin(1);
t1 = t6 * t18;
t2 = mrSges(5,3) * t58;
t3 = mrSges(5,3) * t57;
t4 = t6 * mrSges(4,1);
t9 = -pkin(3) - t11;
t7 = t9 * t17;
t36 = t1 + t2 + t3 + t38 - t4 + t7 - t56;
t27 = Ifges(5,5) * t44;
t24 = t31 * pkin(2) + pkin(7);
t10 = pkin(7) + t12;
t13 = [-0.2e1 * t56 + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 - 0.2e1 * t4 + 0.2e1 * t7 + 0.2e1 * t37 + 0.2e1 * m(5) * (t10 * t42 + t9 * t6) + 0.2e1 * m(4) * (-t11 * t6 + t12 * t5) + t38; t36 + m(5) * (t25 * t6 + (t57 + t58) * t24) + (m(4) * (t31 * t5 - t34 * t6) + (m(5) * (t61 * t10 + t31 * t9) + m(4) * (-t11 * t31 + t12 * t34) + t40) * qJD(3)) * pkin(2) + t37 + t60; 0.2e1 * t15 + 0.2e1 * t20 + 0.2e1 * t21 + 0.2e1 * t8 + 0.2e1 * (m(5) * (t61 * t24 + t25 * t31) + t40) * t48 + t38; m(5) * (-pkin(3) * t6 + pkin(7) * t42) - t59 + t36; -t59 + (m(5) * (-pkin(3) * t31 + t61 * pkin(7)) + t40) * t48 + t38 + t60; t38 - 0.2e1 * t59; t27 - t39 * t5 + (t18 * t10 - t53) * qJD(4); t27 - t39 * t43 + (t18 * t24 - t53) * qJD(4); t27 + (t18 * pkin(7) - t53) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1), t13(2), t13(4), t13(7); t13(2), t13(3), t13(5), t13(8); t13(4), t13(5), t13(6), t13(9); t13(7), t13(8), t13(9), t13(10);];
Mq = res;

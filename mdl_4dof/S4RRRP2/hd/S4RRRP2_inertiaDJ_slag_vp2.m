% Calculate time derivative of joint inertia matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:54
% EndTime: 2019-12-31 17:12:55
% DurationCPUTime: 0.41s
% Computational Cost: add. (253->92), mult. (633->124), div. (0->0), fcn. (335->4), ass. (0->50)
t34 = sin(qJ(3));
t67 = -0.2e1 * t34;
t78 = Ifges(5,4) + Ifges(4,4);
t77 = Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,2);
t36 = cos(qJ(3));
t75 = t34 ^ 2 + t36 ^ 2;
t66 = 0.2e1 * t36;
t37 = cos(qJ(2));
t72 = t75 * t37;
t35 = sin(qJ(2));
t71 = (t75 * mrSges(4,3) - mrSges(3,2)) * t37 + (-t36 * mrSges(4,1) + t34 * mrSges(4,2) - mrSges(3,1)) * t35;
t48 = qJD(3) * t36;
t49 = qJD(3) * t34;
t38 = (t77 * t34 + t78 * t66) * t48 + (t77 * t36 + t78 * t67) * t49;
t70 = 2 * m(5);
t13 = mrSges(5,1) * t49 + mrSges(5,2) * t48;
t69 = 0.2e1 * t13;
t18 = -t36 * mrSges(5,1) + t34 * mrSges(5,2);
t68 = 0.2e1 * t18;
t65 = mrSges(5,3) * pkin(3);
t64 = t37 * pkin(1);
t54 = -Ifges(4,6) - Ifges(5,6);
t53 = -qJ(4) - pkin(6);
t52 = (Ifges(4,5) + Ifges(5,5)) * t48;
t51 = pkin(1) * qJD(2);
t23 = t35 * pkin(1) + pkin(6);
t50 = -qJ(4) - t23;
t47 = 0.2e1 * qJD(3);
t46 = t37 * t51;
t45 = pkin(3) * t49;
t44 = m(5) * pkin(3) + mrSges(5,1);
t25 = -t36 * pkin(3) - pkin(2);
t43 = qJD(3) * t53;
t42 = qJD(3) * t50;
t39 = mrSges(4,1) * t34 + mrSges(4,2) * t36;
t31 = t36 * qJ(4);
t30 = t36 * qJD(4);
t24 = -pkin(2) - t64;
t20 = t36 * pkin(6) + t31;
t17 = t53 * t34;
t16 = t25 - t64;
t15 = t35 * t51 + t45;
t14 = t39 * qJD(3);
t10 = t36 * t23 + t31;
t9 = t50 * t34;
t4 = -t34 * qJD(4) + t36 * t43;
t3 = t34 * t43 + t30;
t2 = (-qJD(4) - t46) * t34 + t36 * t42;
t1 = t34 * t42 + t36 * t46 + t30;
t5 = [(t10 * t1 + t16 * t15 + t9 * t2) * t70 + t15 * t68 + t16 * t69 + 0.2e1 * t24 * t14 + (t1 * t66 + t2 * t67 + (-t10 * t34 - t36 * t9) * t47) * mrSges(5,3) + 0.2e1 * (m(4) * (t23 * t72 + t24 * t35) + t71) * t51 + t38; m(5) * (t20 * t1 + t3 * t10 + t25 * t15 + t16 * t45 + t17 * t2 + t4 * t9) + (t15 + t45) * t18 + (t24 - pkin(2)) * t14 + (t25 + t16) * t13 + (m(4) * (-pkin(2) * t35 + pkin(6) * t72) + t71) * t51 + ((t1 + t3) * t36 + (-t2 - t4) * t34 + ((-t17 - t9) * t36 + (-t10 - t20) * t34) * qJD(3)) * mrSges(5,3) + t38; -0.2e1 * pkin(2) * t14 + t45 * t68 + t25 * t69 + (t17 * t4 + t20 * t3 + t25 * t45) * t70 + (t3 * t66 + t4 * t67 + (-t17 * t36 - t20 * t34) * t47) * mrSges(5,3) + t38; -t1 * mrSges(5,2) + t44 * t2 - t39 * t46 + ((-mrSges(4,1) * t23 - t65) * t36 + (mrSges(4,2) * t23 + t54) * t34) * qJD(3) + t52; -t3 * mrSges(5,2) + t44 * t4 + ((-mrSges(4,1) * pkin(6) - t65) * t36 + (mrSges(4,2) * pkin(6) + t54) * t34) * qJD(3) + t52; 0; m(5) * t15 + t13; m(5) * t45 + t13; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1), t5(2), t5(4), t5(7); t5(2), t5(3), t5(5), t5(8); t5(4), t5(5), t5(6), t5(9); t5(7), t5(8), t5(9), t5(10);];
Mq = res;

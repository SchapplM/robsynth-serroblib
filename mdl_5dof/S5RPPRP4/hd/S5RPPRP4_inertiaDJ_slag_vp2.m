% Calculate time derivative of joint inertia matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:04
% EndTime: 2019-12-31 17:52:05
% DurationCPUTime: 0.31s
% Computational Cost: add. (208->65), mult. (405->102), div. (0->0), fcn. (244->4), ass. (0->36)
t54 = Ifges(5,4) + Ifges(6,4);
t16 = sin(pkin(7));
t30 = t16 * qJD(2);
t18 = sin(qJ(4));
t32 = qJD(4) * t18;
t10 = -pkin(4) * t32 + t30;
t19 = cos(qJ(4));
t8 = (-t18 * mrSges(6,1) - t19 * mrSges(6,2)) * qJD(4);
t51 = m(6) * t10 + t8;
t34 = t18 ^ 2 + t19 ^ 2;
t17 = cos(pkin(7));
t20 = -pkin(1) - pkin(2);
t35 = t17 * qJ(2) + t16 * t20;
t7 = -pkin(6) + t35;
t36 = qJ(5) - t7;
t3 = t36 * t18;
t4 = t36 * t19;
t50 = -t18 * t4 + t19 * t3;
t49 = Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,2);
t48 = 2 * mrSges(6,3);
t33 = qJD(2) * t17;
t23 = -qJD(5) + t33;
t25 = qJD(4) * t36;
t1 = t18 * t25 + t19 * t23;
t45 = t1 * t19;
t2 = -t18 * t23 + t19 * t25;
t44 = t18 * t2;
t37 = mrSges(5,2) + mrSges(6,2);
t27 = m(6) * pkin(4) + mrSges(6,1);
t26 = -t16 * qJ(2) + t17 * t20;
t22 = -mrSges(5,1) - t27;
t6 = pkin(3) - t26;
t21 = -mrSges(5,1) * t18 - mrSges(5,2) * t19;
t9 = t21 * qJD(4);
t5 = t19 * pkin(4) + t6;
t11 = [-0.2e1 * mrSges(6,3) * t45 + 0.2e1 * m(6) * (-t4 * t1 + t5 * t10 + t3 * t2) + 0.2e1 * t10 * (t19 * mrSges(6,1) - t18 * mrSges(6,2)) + 0.2e1 * t5 * t8 + 0.2e1 * t6 * t9 + t44 * t48 + 0.2e1 * (t19 * mrSges(5,1) - t18 * mrSges(5,2) + mrSges(4,1)) * t30 + (-0.2e1 * t54 * t18 + t49 * t19) * t32 + 0.2e1 * (m(5) * (t17 * t34 * t7 + t16 * t6) + m(4) * (-t16 * t26 + t17 * t35) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2) + 0.2e1 * (-t34 * mrSges(5,3) + mrSges(4,2)) * t33 + (t50 * t48 + (t49 * t18 + 0.2e1 * t54 * t19) * t19) * qJD(4); (-t9 - t51) * t17 + (m(6) * (-t50 * qJD(4) - t44 + t45) + m(5) * (-0.1e1 + t34) * t33) * t16; 0; m(6) * (t1 * t18 + t2 * t19 + (-t18 * t3 - t19 * t4) * qJD(4)); 0; 0; -t1 * mrSges(6,2) + t27 * t2 + t21 * t33 + ((mrSges(5,2) * t7 + Ifges(5,6) + Ifges(6,6)) * t18 + (-mrSges(5,1) * t7 + mrSges(6,3) * pkin(4) - Ifges(5,5) - Ifges(6,5)) * t19) * qJD(4); (t18 * t37 + t19 * t22) * t16 * qJD(4); (t18 * t22 - t19 * t37) * qJD(4); 0; t51; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;

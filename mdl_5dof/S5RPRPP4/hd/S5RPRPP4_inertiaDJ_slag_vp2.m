% Calculate time derivative of joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:20
% DurationCPUTime: 0.64s
% Computational Cost: add. (455->89), mult. (961->130), div. (0->0), fcn. (711->4), ass. (0->45)
t77 = mrSges(5,3) + mrSges(6,2);
t29 = sin(pkin(7));
t53 = cos(pkin(7));
t60 = cos(qJ(3));
t36 = t53 * t60;
t30 = sin(qJ(3));
t52 = qJD(3) * t30;
t15 = -qJD(3) * t36 + t29 * t52;
t19 = t29 * t60 + t30 * t53;
t16 = t19 * qJD(3);
t18 = t29 * t30 - t36;
t31 = -pkin(1) - pkin(6);
t54 = qJ(4) - t31;
t33 = t54 * t60;
t12 = -qJD(3) * t33 - t30 * qJD(4);
t32 = -qJD(4) * t60 + t52 * t54;
t5 = t29 * t12 - t32 * t53;
t6 = t12 * t53 + t29 * t32;
t20 = t54 * t30;
t8 = -t29 * t20 + t33 * t53;
t9 = -t20 * t53 - t29 * t33;
t76 = -t15 * t9 + t16 * t8 + t18 * t5 + t19 * t6;
t75 = m(5) + m(6);
t57 = t18 * t16;
t10 = t19 * t15;
t74 = -mrSges(5,1) - mrSges(6,1);
t73 = -Ifges(4,1) + Ifges(4,2);
t72 = t10 - t57;
t39 = qJD(3) * t60;
t71 = -mrSges(4,1) * t52 - mrSges(4,2) * t39;
t70 = t15 * t29 + t53 * t16;
t23 = t29 * pkin(3) + qJ(5);
t25 = -pkin(3) * t53 - pkin(4);
t69 = qJD(5) * t19 - t23 * t15 + t25 * t16;
t68 = m(6) * t23 + mrSges(6,3);
t65 = 2 * qJD(2);
t64 = m(5) * pkin(3);
t26 = t30 * pkin(3) + qJ(2);
t21 = pkin(3) * t39 + qJD(2);
t45 = t8 * t5 + t9 * t6;
t42 = -t15 * mrSges(5,1) - t16 * mrSges(5,2);
t41 = -t15 * mrSges(6,1) + t16 * mrSges(6,3);
t7 = t19 * pkin(4) + t18 * qJ(5) + t26;
t3 = -t15 * pkin(4) + t16 * qJ(5) + t18 * qJD(5) + t21;
t1 = [0.2e1 * m(5) * (t26 * t21 + t45) + 0.2e1 * m(6) * (t7 * t3 + t45) + 0.2e1 * t3 * (t19 * mrSges(6,1) + t18 * mrSges(6,3)) + 0.2e1 * t21 * (t19 * mrSges(5,1) - t18 * mrSges(5,2)) + 0.2e1 * t7 * t41 + 0.2e1 * t26 * t42 + (0.2e1 * Ifges(4,4) * t30 + t73 * t60) * t52 + (-0.2e1 * Ifges(4,4) * t60 + t73 * t30) * t39 + (t30 * mrSges(4,1) + mrSges(4,2) * t60 + mrSges(3,3)) * t65 - 0.2e1 * t77 * t76 + (0.2e1 * (mrSges(4,1) * t60 - mrSges(4,2) * t30) * qJD(3) + ((m(4) + m(3)) * t65)) * qJ(2) + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t57 + 0.2e1 * (-Ifges(5,2) - Ifges(6,3)) * t10 + 0.2e1 * (Ifges(6,5) - Ifges(5,4)) * (t15 * t18 - t16 * t19); 0.2e1 * t77 * t72 + t75 * t76; -0.2e1 * t75 * t72; m(6) * qJD(5) * t9 - Ifges(4,5) * t52 - Ifges(4,6) * t39 + t70 * mrSges(5,3) * pkin(3) + (t29 * t64 - mrSges(5,2) + t68) * t6 + (m(6) * t25 - t53 * t64 + t74) * t5 + t71 * t31 + (-Ifges(6,4) - Ifges(5,5)) * t16 + (-Ifges(6,6) + Ifges(5,6)) * t15 - t69 * mrSges(6,2); -t70 * t64 + m(6) * t69 + t74 * t16 + (mrSges(5,2) - mrSges(6,3)) * t15 + t71; 0.2e1 * t68 * qJD(5); m(5) * t21 + m(6) * t3 + t41 + t42; 0; 0; 0; m(6) * t5 - t16 * mrSges(6,2); m(6) * t16; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

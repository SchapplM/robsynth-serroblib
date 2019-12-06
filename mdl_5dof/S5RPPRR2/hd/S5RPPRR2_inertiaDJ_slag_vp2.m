% Calculate time derivative of joint inertia matrix for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:37
% DurationCPUTime: 0.62s
% Computational Cost: add. (1071->98), mult. (2099->162), div. (0->0), fcn. (1984->6), ass. (0->51)
t50 = sin(pkin(8));
t51 = cos(pkin(8));
t76 = sin(qJ(4));
t63 = qJD(4) * t76;
t55 = cos(qJ(4));
t70 = qJD(4) * t55;
t32 = -t50 * t63 + t51 * t70;
t38 = -t55 * t50 - t76 * t51;
t33 = t38 * qJD(4);
t53 = sin(qJ(5));
t54 = cos(qJ(5));
t85 = t76 * t50 - t55 * t51;
t79 = t38 * t53 - t54 * t85;
t10 = qJD(5) * t79 + t32 * t54 + t33 * t53;
t21 = t38 * t54 + t53 * t85;
t90 = t10 * t21;
t52 = -pkin(1) - qJ(3);
t77 = -pkin(6) + t52;
t40 = t77 * t50;
t41 = t77 * t51;
t16 = t38 * qJD(3) - t40 * t63 + t41 * t70;
t14 = -t32 * pkin(7) + t16;
t25 = t55 * t40 + t76 * t41;
t17 = t85 * qJD(3) - t25 * qJD(4);
t15 = -t33 * pkin(7) + t17;
t24 = -t76 * t40 + t55 * t41;
t18 = pkin(7) * t85 + t24;
t19 = pkin(7) * t38 + t25;
t4 = t18 * t54 - t19 * t53;
t2 = t4 * qJD(5) + t14 * t54 + t15 * t53;
t5 = t18 * t53 + t19 * t54;
t3 = -t5 * qJD(5) - t14 * t53 + t15 * t54;
t56 = t21 * qJD(5) - t32 * t53 + t54 * t33;
t89 = t10 * t5 - t2 * t21 + t3 * t79 + t4 * t56;
t88 = (t21 * t54 + t53 * t79) * qJD(5) - t10 * t53 - t54 * t56;
t80 = t56 * t79;
t62 = (t50 ^ 2 + t51 ^ 2) * qJD(3);
t75 = t32 * t38;
t74 = t33 * t85;
t44 = t50 * pkin(3) + qJ(2);
t69 = 2 * mrSges(5,3);
t67 = mrSges(6,1) * t10 + mrSges(6,2) * t56;
t66 = mrSges(6,1) * t56 - t10 * mrSges(6,2);
t65 = t32 * mrSges(5,1) + t33 * mrSges(5,2);
t61 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t56 - Ifges(6,6) * t10;
t60 = t74 + t75;
t57 = t16 * t38 + t17 * t85 - t24 * t33 - t25 * t32;
t37 = (-mrSges(6,1) * t53 - mrSges(6,2) * t54) * qJD(5) * pkin(4);
t27 = pkin(4) * t32 + qJD(2);
t26 = -pkin(4) * t38 + t44;
t1 = [-0.2e1 * Ifges(5,2) * t75 - 0.2e1 * Ifges(5,1) * t74 + 0.2e1 * t44 * t65 - 0.2e1 * Ifges(6,2) * t90 + 0.2e1 * Ifges(6,1) * t80 + 0.2e1 * t26 * t67 + 0.2e1 * t27 * (-t21 * mrSges(6,1) + mrSges(6,2) * t79) + 0.2e1 * m(5) * (qJD(2) * t44 + t16 * t25 + t17 * t24) + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t52 * t62) + 0.2e1 * m(6) * (t2 * t5 + t26 * t27 + t3 * t4) + t57 * t69 + 0.2e1 * (m(3) * qJ(2) + t50 * mrSges(4,1) - t38 * mrSges(5,1) + t51 * mrSges(4,2) - mrSges(5,2) * t85 + mrSges(3,3)) * qJD(2) + 0.2e1 * mrSges(4,3) * t62 + 0.2e1 * (-t10 * t79 + t21 * t56) * Ifges(6,4) + 0.2e1 * (t32 * t85 + t33 * t38) * Ifges(5,4) - 0.2e1 * t89 * mrSges(6,3); t60 * t69 + m(6) * t89 - m(5) * t57 - m(4) * t62 + (-0.2e1 * t80 + 0.2e1 * t90) * mrSges(6,3); -0.2e1 * m(5) * t60 + 0.2e1 * m(6) * (t80 - t90); m(6) * t27 + (m(5) + m(4)) * qJD(2) + t65 + t67; 0; 0; t17 * mrSges(5,1) - t16 * mrSges(5,2) + Ifges(5,5) * t33 - Ifges(5,6) * t32 + (m(6) * (t2 * t53 + t3 * t54 + (-t4 * t53 + t5 * t54) * qJD(5)) + t88 * mrSges(6,3)) * pkin(4) + t61; -m(6) * t88 * pkin(4) + t33 * mrSges(5,1) - t32 * mrSges(5,2) + t66; 0; 0.2e1 * t37; t61; t66; 0; t37; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

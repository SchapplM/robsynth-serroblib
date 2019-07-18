% Calculate time derivative of joint inertia matrix for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:12
% EndTime: 2019-07-18 13:30:13
% DurationCPUTime: 0.29s
% Computational Cost: add. (288->69), mult. (829->110), div. (0->0), fcn. (474->6), ass. (0->53)
t31 = sin(qJ(4));
t65 = t31 * pkin(3);
t30 = sin(qJ(5));
t28 = t30 ^ 2;
t33 = cos(qJ(5));
t29 = t33 ^ 2;
t50 = t28 + t29;
t35 = cos(qJ(3));
t24 = t35 * pkin(2) + pkin(3);
t32 = sin(qJ(3));
t34 = cos(qJ(4));
t53 = t32 * t34;
t10 = pkin(2) * t53 + t31 * t24;
t8 = pkin(6) + t10;
t64 = t50 * t8;
t63 = Ifges(6,1) - Ifges(6,2);
t52 = t33 * mrSges(6,1);
t16 = t30 * mrSges(6,2) - t52;
t49 = qJD(4) * t31;
t13 = pkin(3) * t16 * t49;
t48 = qJD(4) * t34;
t45 = pkin(3) * t48;
t42 = mrSges(6,3) * t45;
t18 = t28 * t42;
t19 = t29 * t42;
t62 = t13 + t18 + t19;
t6 = t24 * t49 + (t32 * t48 + (t31 * t35 + t53) * qJD(3)) * pkin(2);
t54 = t31 * t32;
t9 = pkin(2) * t54 - t34 * t24;
t61 = t9 * t6;
t5 = t24 * t48 + (-t32 * t49 + (t34 * t35 - t54) * qJD(3)) * pkin(2);
t60 = mrSges(6,3) * t5;
t56 = Ifges(6,6) * t30;
t55 = t31 * mrSges(5,1);
t46 = qJD(5) * t33;
t47 = qJD(5) * t30;
t15 = -mrSges(6,1) * t47 - mrSges(6,2) * t46;
t51 = t34 * t15;
t43 = m(6) * t50;
t41 = -mrSges(6,1) * t30 - mrSges(6,2) * t33;
t40 = pkin(6) * t43 - mrSges(5,2);
t39 = -t34 * t6 + t9 * t49;
t38 = (-0.2e1 * Ifges(6,4) * t30 + t63 * t33) * t47 + (0.2e1 * Ifges(6,4) * t33 + t63 * t30) * t46;
t37 = (-mrSges(4,1) * t32 - mrSges(4,2) * t35) * qJD(3) * pkin(2);
t1 = t6 * t16;
t2 = t28 * t60;
t3 = t29 * t60;
t4 = t6 * mrSges(5,1);
t7 = t9 * t15;
t36 = t1 + t2 + t3 + t38 - t4 - t7;
t27 = Ifges(6,5) * t46;
t23 = pkin(6) + t65;
t11 = [0; 0; -0.2e1 * t5 * mrSges(5,2) + 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t3 - 0.2e1 * t4 - 0.2e1 * t7 + 0.2e1 * t37 + 0.2e1 * m(6) * (t5 * t64 + t61) + 0.2e1 * m(5) * (t10 * t5 + t61) + t38; 0; t36 + (t23 * t43 - mrSges(5,2)) * t5 + (t51 + (-t34 * mrSges(5,2) - t55) * qJD(4) + m(6) * (t48 * t64 + t39) + m(5) * (t10 * t48 + t31 * t5 + t39)) * pkin(3) + t37 + t62; 0.2e1 * t13 + 0.2e1 * t18 + 0.2e1 * t19 + t38 + 0.2e1 * (t51 + (-t55 + (-mrSges(5,2) + (t23 * t50 - t65) * m(6)) * t34) * qJD(4)) * pkin(3); 0; t40 * t5 + t36; (t51 + (t40 * t34 - t55) * qJD(4)) * pkin(3) + t38 + t62; t38; t15; t27 + t41 * t5 + (-t8 * t52 + (mrSges(6,2) * t8 - Ifges(6,6)) * t30) * qJD(5); t27 + t41 * t45 + (t16 * t23 - t56) * qJD(5); t27 + (t16 * pkin(6) - t56) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq  = res;

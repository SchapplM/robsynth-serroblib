% Calculate time derivative of joint inertia matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:04
% EndTime: 2019-12-05 17:47:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (1152->116), mult. (2273->189), div. (0->0), fcn. (2042->6), ass. (0->63)
t57 = sin(pkin(8));
t61 = cos(qJ(3));
t81 = cos(pkin(8));
t70 = t81 * t61;
t59 = sin(qJ(3));
t80 = qJD(3) * t59;
t40 = -qJD(3) * t70 + t57 * t80;
t45 = -t57 * t61 - t81 * t59;
t41 = t45 * qJD(3);
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t46 = -t57 * t59 + t70;
t95 = t45 * t58 + t60 * t46;
t10 = qJD(5) * t95 - t40 * t60 + t58 * t41;
t21 = t45 * t60 - t58 * t46;
t103 = t10 * t21;
t62 = -pkin(1) - pkin(6);
t82 = qJ(4) - t62;
t34 = -t61 * qJD(4) + t82 * t80;
t48 = t82 * t61;
t35 = -qJD(3) * t48 - t59 * qJD(4);
t16 = t81 * t34 - t35 * t57;
t14 = -pkin(7) * t41 + t16;
t17 = t57 * t34 + t81 * t35;
t15 = pkin(7) * t40 + t17;
t47 = t82 * t59;
t24 = t47 * t57 - t81 * t48;
t18 = -pkin(7) * t46 + t24;
t25 = -t81 * t47 - t57 * t48;
t19 = pkin(7) * t45 + t25;
t4 = t18 * t60 - t58 * t19;
t2 = t4 * qJD(5) + t58 * t14 + t15 * t60;
t5 = t58 * t18 + t19 * t60;
t3 = -t5 * qJD(5) + t14 * t60 - t58 * t15;
t63 = t21 * qJD(5) + t58 * t40 + t60 * t41;
t102 = t10 * t5 - t2 * t21 + t3 * t95 + t4 * t63;
t53 = t81 * pkin(3) + pkin(4);
t88 = pkin(3) * t57;
t37 = t53 * t60 - t58 * t88;
t31 = t37 * qJD(5);
t38 = t58 * t53 + t60 * t88;
t32 = t38 * qJD(5);
t101 = t10 * t38 - t31 * t21 - t32 * t95 + t37 * t63;
t79 = qJD(3) * t61;
t96 = t63 * t95;
t94 = -mrSges(4,1) * t80 - mrSges(4,2) * t79;
t93 = t40 * t57 - t81 * t41;
t90 = m(5) * pkin(3);
t86 = t41 * t46;
t85 = t45 * t40;
t54 = t59 * pkin(3) + qJ(2);
t49 = pkin(3) * t79 + qJD(2);
t78 = 2 * mrSges(5,3);
t75 = mrSges(6,1) * t10 + mrSges(6,2) * t63;
t74 = mrSges(6,1) * t63 - t10 * mrSges(6,2);
t73 = -t40 * mrSges(5,1) + t41 * mrSges(5,2);
t69 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t63 - Ifges(6,6) * t10;
t68 = -t32 * mrSges(6,1) - t31 * mrSges(6,2);
t67 = t85 + t86;
t64 = t16 * t46 - t17 * t45 + t24 * t41 - t25 * t40;
t30 = -pkin(4) * t45 + t54;
t26 = -pkin(4) * t40 + t49;
t1 = [0.2e1 * t54 * t73 + 0.2e1 * Ifges(5,2) * t85 + 0.2e1 * Ifges(5,1) * t86 + 0.2e1 * t49 * (-t45 * mrSges(5,1) + t46 * mrSges(5,2)) + 0.2e1 * Ifges(6,1) * t96 - 0.2e1 * Ifges(6,2) * t103 + 0.2e1 * t26 * (-t21 * mrSges(6,1) + mrSges(6,2) * t95) + 0.2e1 * t30 * t75 + 0.2e1 * m(5) * (t16 * t24 + t17 * t25 + t49 * t54) + 0.2e1 * m(6) * (t2 * t5 + t26 * t30 + t3 * t4) + 0.2e1 * (-t10 * t95 + t21 * t63) * Ifges(6,4) + 0.2e1 * (t40 * t46 + t45 * t41) * Ifges(5,4) - 0.2e1 * t102 * mrSges(6,3) - t64 * t78 + 0.2e1 * (qJ(2) * (mrSges(4,1) * t61 - mrSges(4,2) * t59) + (t59 ^ 2 - t61 ^ 2) * Ifges(4,4)) * qJD(3) + 0.2e1 * (t59 * mrSges(4,1) + mrSges(4,2) * t61 + mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * qJD(2) + 0.2e1 * (-Ifges(4,1) + Ifges(4,2)) * t59 * t79; -t67 * t78 + m(6) * t102 + m(5) * t64 + (-0.2e1 * t96 + 0.2e1 * t103) * mrSges(6,3); 0.2e1 * m(5) * t67 + 0.2e1 * m(6) * (t96 - t103); m(6) * (t2 * t38 + t3 * t37 + t31 * t5 - t32 * t4) + t69 + Ifges(5,6) * t40 + Ifges(5,5) * t41 + t16 * mrSges(5,1) - t17 * mrSges(5,2) - Ifges(4,5) * t80 - Ifges(4,6) * t79 + (t81 * t16 + t17 * t57) * t90 + t93 * mrSges(5,3) * pkin(3) + t94 * t62 - t101 * mrSges(6,3); m(6) * t101 + t41 * mrSges(5,1) + t40 * mrSges(5,2) - t93 * t90 + t74 + t94; 0.2e1 * m(6) * (t31 * t38 - t32 * t37) + 0.2e1 * t68; m(5) * t49 + m(6) * t26 + t73 + t75; 0; 0; 0; t69; t74; t68; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

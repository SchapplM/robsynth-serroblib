% Calculate time derivative of joint inertia matrix for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:38
% EndTime: 2019-12-05 16:45:40
% DurationCPUTime: 0.58s
% Computational Cost: add. (440->99), mult. (1192->144), div. (0->0), fcn. (808->6), ass. (0->60)
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t82 = t52 ^ 2 + t55 ^ 2;
t104 = mrSges(6,2) + mrSges(5,3);
t103 = 2 * Ifges(5,4) - 2 * Ifges(6,5);
t36 = -t55 * mrSges(6,1) - t52 * mrSges(6,3);
t37 = -t55 * mrSges(5,1) + t52 * mrSges(5,2);
t102 = t36 + t37;
t56 = cos(qJ(3));
t80 = qJD(3) * t56;
t75 = pkin(2) * t80;
t101 = t82 * t75;
t53 = sin(qJ(3));
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t30 = t53 * t54 - t56 * t57;
t98 = qJD(2) + qJD(3);
t13 = t98 * t30;
t100 = t82 * t13;
t99 = Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3);
t78 = qJD(4) * t55;
t97 = (t104 * t82 - mrSges(4,2)) * t75;
t66 = t52 * mrSges(6,1) - t55 * mrSges(6,3);
t32 = t66 * qJD(4);
t96 = 0.2e1 * t32;
t95 = m(5) / 0.2e1;
t94 = t56 * pkin(2);
t31 = t53 * t57 + t56 * t54;
t14 = t98 * t31;
t87 = t30 * t14;
t86 = t31 * t13;
t85 = t100 * pkin(7);
t42 = t53 * pkin(2) + pkin(7);
t84 = t101 * t42;
t83 = t101 * pkin(7);
t81 = qJD(3) * t53;
t79 = qJD(4) * t52;
t77 = qJD(5) * t55;
t76 = m(6) * t77;
t74 = pkin(2) * t81;
t73 = t30 * t81;
t72 = -t100 * t42 + t101 * t31;
t71 = mrSges(6,2) * t77 + Ifges(6,6) * t79 + (Ifges(6,4) + Ifges(5,5)) * t78;
t67 = t52 * mrSges(5,1) + t55 * mrSges(5,2);
t65 = -t55 * pkin(4) - t52 * qJ(5);
t64 = (-mrSges(4,1) + t37) * t74;
t35 = -pkin(3) + t65;
t20 = pkin(4) * t79 - qJ(5) * t78 - t52 * qJD(5);
t63 = t65 * mrSges(6,2) - Ifges(5,6) * t52;
t33 = t67 * qJD(4);
t62 = t13 * mrSges(4,2) + (t32 + t33) * t30 + (-mrSges(4,1) + t102) * t14 - t104 * t100;
t61 = (-t103 * t52 + t99 * t55) * t79 + (t103 * t55 + t99 * t52) * t78;
t60 = m(6) * t65 + t102;
t59 = m(6) * (-pkin(4) * t52 + qJ(5) * t55) - t66 - t67;
t58 = qJD(4) * t60 + t76;
t45 = mrSges(6,2) * t78;
t43 = -pkin(3) - t94;
t24 = t35 - t94;
t17 = t20 + t74;
t1 = [0.2e1 * m(4) * (-t86 + t87) + 0.4e1 * (t95 + m(6) / 0.2e1) * (-t82 * t86 + t87); t62 + m(5) * (t43 * t14 + t72) + 0.2e1 * (m(4) * (-t13 * t53 - t14 * t56 + t31 * t80 + t73) / 0.2e1 + t73 * t95) * pkin(2) + m(6) * (t24 * t14 + t17 * t30 + t72) + (-t54 * mrSges(3,1) - t57 * mrSges(3,2)) * qJD(2); 0.2e1 * t17 * t36 + t24 * t96 + 0.2e1 * t43 * t33 + 0.2e1 * t64 + 0.2e1 * m(6) * (t24 * t17 + t84) + 0.2e1 * m(5) * (t43 * t74 + t84) + 0.2e1 * t97 + t61; m(5) * (-pkin(3) * t14 - t85) + m(6) * (t35 * t14 + t20 * t30 - t85) + t62; (t17 + t20) * t36 + (t43 - pkin(3)) * t33 + (t24 + t35) * t32 + t64 + m(6) * (t35 * t17 + t20 * t24 + t83) + m(5) * (-pkin(3) * t74 + t83) + t97 + t61; -0.2e1 * pkin(3) * t33 + t35 * t96 + 0.2e1 * (m(6) * t35 + t36) * t20 + t61; -t13 * t59 + t31 * t58; t42 * t76 + t59 * t75 + (t42 * t60 + t63) * qJD(4) + t71; pkin(7) * t58 + qJD(4) * t63 + t71; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); (-t52 * t13 + t31 * t78) * m(6); t45 + (t42 * t78 + t52 * t75) * m(6); m(6) * pkin(7) * t78 + t45; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

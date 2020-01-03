% Calculate time derivative of joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:53
% EndTime: 2019-12-31 20:54:56
% DurationCPUTime: 1.06s
% Computational Cost: add. (1666->164), mult. (3882->252), div. (0->0), fcn. (3307->6), ass. (0->77)
t84 = (mrSges(6,2) + mrSges(5,3));
t102 = 2 * t84;
t101 = -mrSges(5,1) - mrSges(6,1);
t63 = sin(qJ(3));
t64 = sin(qJ(2));
t65 = cos(qJ(3));
t66 = cos(qJ(2));
t48 = t63 * t66 + t65 * t64;
t94 = -pkin(7) - pkin(6);
t52 = t94 * t64;
t49 = t63 * t52;
t53 = t94 * t66;
t34 = -t65 * t53 + t49;
t100 = qJD(2) + qJD(3);
t99 = 2 * m(5);
t98 = 2 * m(6);
t85 = t65 * t66;
t47 = -t63 * t64 + t85;
t59 = -pkin(2) * t66 - pkin(1);
t35 = -t47 * pkin(3) + t59;
t97 = 0.2e1 * t35;
t96 = 0.2e1 * t59;
t95 = m(5) * pkin(3);
t62 = sin(pkin(8));
t93 = pkin(3) * t62;
t24 = qJ(4) * t47 + t34;
t33 = t52 * t65 + t53 * t63;
t69 = -t48 * qJ(4) + t33;
t82 = cos(pkin(8));
t10 = t24 * t62 - t82 * t69;
t73 = t82 * t63;
t83 = pkin(2) * qJD(3);
t39 = (t62 * t65 + t73) * t83;
t92 = t10 * t39;
t31 = t100 * t47;
t32 = t100 * t48;
t19 = t82 * t31 - t62 * t32;
t91 = t19 * mrSges(6,2);
t28 = t62 * t47 + t82 * t48;
t90 = t28 * t39;
t80 = qJD(3) * t65;
t77 = pkin(2) * t80;
t81 = qJD(3) * t63;
t78 = pkin(2) * t81;
t40 = -t62 * t78 + t82 * t77;
t36 = qJD(5) + t40;
t89 = t36 * mrSges(6,3);
t88 = t40 * mrSges(5,2);
t58 = pkin(2) * t65 + pkin(3);
t42 = pkin(2) * t73 + t62 * t58;
t79 = 0.2e1 * t66;
t11 = t82 * t24 + t62 * t69;
t68 = (t94 * t85 - t49) * qJD(2);
t22 = -qJD(3) * t34 + t68;
t71 = -t31 * qJ(4) - t48 * qJD(4);
t21 = t48 * qJD(2) * t94 + t52 * t80 + t53 * t81;
t8 = -t32 * qJ(4) + t47 * qJD(4) + t21;
t4 = t62 * t8 - t82 * (t22 + t71);
t5 = t82 * t8 + (-t52 * t81 + t53 * t80 + t68 + t71) * t62;
t76 = t10 * t4 + t11 * t5;
t25 = qJD(2) * t64 * pkin(2) + t32 * pkin(3);
t74 = t82 * pkin(3);
t72 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t41 = -t62 * t63 * pkin(2) + t82 * t58;
t18 = t31 * t62 + t82 * t32;
t67 = t22 * mrSges(4,1) - t21 * mrSges(4,2) + Ifges(4,5) * t31 - Ifges(4,6) * t32 + (-mrSges(5,2) + mrSges(6,3)) * t5 + t101 * t4 + (Ifges(6,4) + Ifges(5,5)) * t19 + (-Ifges(5,6) + Ifges(6,6)) * t18;
t61 = qJD(5) * mrSges(6,3);
t57 = -t74 - pkin(4);
t56 = qJ(5) + t93;
t38 = -pkin(4) - t41;
t37 = qJ(5) + t42;
t27 = -t82 * t47 + t48 * t62;
t13 = t19 * mrSges(5,2);
t12 = t18 * mrSges(6,1);
t9 = t27 * pkin(4) - t28 * qJ(5) + t35;
t6 = t18 * pkin(4) - t19 * qJ(5) - t28 * qJD(5) + t25;
t1 = [(mrSges(4,1) * t32 + mrSges(4,2) * t31) * t96 - 0.2e1 * t47 * Ifges(4,2) * t32 + 0.2e1 * t31 * t48 * Ifges(4,1) + 0.2e1 * t6 * (mrSges(6,1) * t27 - mrSges(6,3) * t28) + 0.2e1 * t25 * (mrSges(5,1) * t27 + mrSges(5,2) * t28) + 0.2e1 * t9 * t12 + t13 * t97 + (t25 * t35 + t76) * t99 + (t6 * t9 + t76) * t98 + 0.2e1 * m(4) * (t34 * t21 + t33 * t22) + 0.2e1 * (t47 * t31 - t32 * t48) * Ifges(4,4) + 0.2e1 * (t21 * t47 - t22 * t48 - t33 * t31 - t34 * t32) * mrSges(4,3) + (-0.2e1 * t9 * mrSges(6,3) + t27 * t72 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t28 + t10 * t102) * t19 + (mrSges(5,1) * t97 + t28 * t72 + 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t27 - 0.2e1 * t84 * t11) * t18 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t66) * t79 + (-0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t96 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t47 + mrSges(4,2) * t48) - 0.2e1 * Ifges(3,4) * t64 + (Ifges(3,1) - Ifges(3,2)) * t79) * t64) * qJD(2) + (-t5 * t27 + t4 * t28) * t102; m(5) * (t40 * t11 - t4 * t41 + t42 * t5 + t92) + m(6) * (t11 * t36 + t37 * t5 + t38 * t4 + t92) + (-t18 * t42 - t19 * t41 - t27 * t40 + t90) * mrSges(5,3) + (-t18 * t37 + t19 * t38 - t27 * t36 + t90) * mrSges(6,2) + (m(4) * (t21 * t63 + t22 * t65 + (-t33 * t63 + t34 * t65) * qJD(3)) + (-t65 * t31 - t63 * t32 + (t47 * t65 + t48 * t63) * qJD(3)) * mrSges(4,3)) * pkin(2) + (Ifges(3,5) * t66 - Ifges(3,6) * t64 + (-mrSges(3,1) * t66 + mrSges(3,2) * t64) * pkin(6)) * qJD(2) + t67; -0.2e1 * t88 + 0.2e1 * t89 + (-t39 * t41 + t40 * t42) * t99 + (t36 * t37 + t38 * t39) * t98 + 0.2e1 * t101 * t39 + 0.2e1 * (-mrSges(4,1) * t63 - mrSges(4,2) * t65) * t83; (-t82 * t4 + t5 * t62) * t95 + m(6) * (qJD(5) * t11 + t4 * t57 + t5 * t56) + t57 * t91 + t67 + (-t18 * t93 - t19 * t74) * mrSges(5,3) + (-qJD(5) * t27 - t56 * t18) * mrSges(6,2); -mrSges(4,2) * t77 - mrSges(4,1) * t78 + t40 * t62 * t95 + t61 + m(6) * (qJD(5) * t37 + t36 * t56) + t89 - t88 + (m(6) * t57 - t82 * t95 + t101) * t39; 0.2e1 * m(6) * qJD(5) * t56 + 0.2e1 * t61; m(5) * t25 + m(6) * t6 + t18 * mrSges(5,1) - t19 * mrSges(6,3) + t12 + t13; 0; 0; 0; m(6) * t4 + t91; m(6) * t39; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

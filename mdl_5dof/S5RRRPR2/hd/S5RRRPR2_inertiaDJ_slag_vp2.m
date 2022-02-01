% Calculate time derivative of joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:34
% DurationCPUTime: 0.39s
% Computational Cost: add. (683->105), mult. (1827->165), div. (0->0), fcn. (1183->8), ass. (0->77)
t90 = Ifges(6,1) - Ifges(6,2);
t57 = cos(qJ(3));
t44 = t57 * pkin(2) + pkin(3);
t52 = cos(pkin(9));
t51 = sin(pkin(9));
t54 = sin(qJ(3));
t77 = t51 * t54;
t29 = -pkin(2) * t77 + t52 * t44;
t24 = -pkin(4) - t29;
t56 = cos(qJ(5));
t67 = qJD(5) * t56;
t53 = sin(qJ(5));
t68 = qJD(5) * t53;
t36 = -mrSges(6,1) * t68 - mrSges(6,2) * t67;
t16 = t24 * t36;
t71 = pkin(2) * qJD(3);
t76 = t52 * t54;
t26 = (t51 * t57 + t76) * t71;
t73 = t56 * mrSges(6,1);
t37 = t53 * mrSges(6,2) - t73;
t17 = t26 * t37;
t27 = (t52 * t57 - t77) * t71;
t49 = t53 ^ 2;
t84 = mrSges(6,3) * t49;
t18 = t27 * t84;
t50 = t56 ^ 2;
t83 = mrSges(6,3) * t50;
t19 = t27 * t83;
t22 = t26 * mrSges(5,1);
t89 = t17 + t18 + t19 - t16 - t22;
t88 = 2 * m(5);
t87 = 2 * m(6);
t86 = m(5) * pkin(3);
t58 = cos(qJ(2));
t45 = t58 * pkin(1) + pkin(2);
t55 = sin(qJ(2));
t69 = qJD(3) * t57;
t70 = qJD(3) * t54;
t75 = t54 * t55;
t14 = t45 * t69 + (-t55 * t70 + (t57 * t58 - t75) * qJD(2)) * pkin(1);
t74 = t55 * t57;
t15 = -t45 * t70 + (-t55 * t69 + (-t54 * t58 - t74) * qJD(2)) * pkin(1);
t6 = t52 * t14 + t51 * t15;
t85 = t6 * mrSges(5,2);
t80 = Ifges(6,6) * t53;
t79 = t14 * mrSges(4,2);
t78 = t27 * mrSges(5,2);
t31 = -pkin(1) * t75 + t57 * t45;
t28 = pkin(3) + t31;
t32 = pkin(1) * t74 + t54 * t45;
t11 = t51 * t28 + t52 * t32;
t30 = pkin(2) * t76 + t51 * t44;
t72 = t49 + t50;
t66 = t72 * t6;
t65 = t72 * t27;
t42 = t51 * pkin(3) + pkin(8);
t64 = t72 * t42;
t63 = -mrSges(6,1) * t53 - mrSges(6,2) * t56;
t10 = t52 * t28 - t51 * t32;
t62 = (-0.2e1 * Ifges(6,4) * t53 + t90 * t56) * t68 + (0.2e1 * Ifges(6,4) * t56 + t90 * t53) * t67;
t61 = (-mrSges(3,1) * t55 - mrSges(3,2) * t58) * qJD(2) * pkin(1);
t60 = (-mrSges(4,1) * t54 - mrSges(4,2) * t57) * t71;
t5 = t51 * t14 - t52 * t15;
t1 = t5 * t37;
t13 = t15 * mrSges(4,1);
t2 = t6 * t84;
t3 = t6 * t83;
t4 = t5 * mrSges(5,1);
t8 = -pkin(4) - t10;
t7 = t8 * t36;
t59 = t1 + t13 + t2 + t3 - t4 + t62 - t7 - t79;
t48 = Ifges(6,5) * t67;
t43 = -t52 * pkin(3) - pkin(4);
t25 = pkin(8) + t30;
t23 = t43 * t36;
t9 = pkin(8) + t11;
t12 = [-0.2e1 * t79 - 0.2e1 * t85 + 0.2e1 * t1 + 0.2e1 * t13 + 0.2e1 * t2 + 0.2e1 * t3 - 0.2e1 * t4 - 0.2e1 * t7 + 0.2e1 * t61 + (-t10 * t5 + t11 * t6) * t88 + (t8 * t5 + t9 * t66) * t87 + 0.2e1 * m(4) * (t32 * t14 + t31 * t15) + t62; t59 + m(6) * (t24 * t5 + t25 * t66 + t26 * t8 + t9 * t65) + (m(4) * (t14 * t54 + t15 * t57 - t31 * t70 + t32 * t69) - mrSges(4,2) * t69 - mrSges(4,1) * t70) * pkin(2) + m(5) * (-t26 * t10 + t27 * t11 - t29 * t5 + t30 * t6) + (-t27 - t6) * mrSges(5,2) + t61 + t89; -0.2e1 * t78 - 0.2e1 * t16 + 0.2e1 * t17 + 0.2e1 * t18 + 0.2e1 * t19 - 0.2e1 * t22 + 0.2e1 * t60 + (t24 * t26 + t25 * t65) * t87 + (-t29 * t26 + t30 * t27) * t88 + t62; t59 + m(6) * (t43 * t5 + t6 * t64) + (-t5 * t52 + t51 * t6) * t86 - t23 - t85; -t78 - t23 + t60 + m(6) * (t43 * t26 + t27 * t64) + (-t26 * t52 + t27 * t51) * t86 + t62 + t89; -0.2e1 * t23 + t62; 0; 0; 0; 0; t48 + t63 * t6 + (-t9 * t73 + (mrSges(6,2) * t9 - Ifges(6,6)) * t53) * qJD(5); t48 + t63 * t27 + (t37 * t25 - t80) * qJD(5); t48 + (t37 * t42 - t80) * qJD(5); t36; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t12(1), t12(2), t12(4), t12(7), t12(11); t12(2), t12(3), t12(5), t12(8), t12(12); t12(4), t12(5), t12(6), t12(9), t12(13); t12(7), t12(8), t12(9), t12(10), t12(14); t12(11), t12(12), t12(13), t12(14), t12(15);];
Mq = res;

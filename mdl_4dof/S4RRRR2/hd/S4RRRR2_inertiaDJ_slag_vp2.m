% Calculate time derivative of joint inertia matrix for
% S4RRRR2
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:10
% DurationCPUTime: 0.54s
% Computational Cost: add. (714->115), mult. (1682->175), div. (0->0), fcn. (1243->6), ass. (0->69)
t51 = sin(qJ(3));
t54 = cos(qJ(3));
t96 = t51 ^ 2 + t54 ^ 2;
t95 = Ifges(4,1) - Ifges(4,2);
t52 = sin(qJ(2));
t43 = t52 * pkin(1) + pkin(6);
t83 = -pkin(7) - t43;
t30 = t83 * t51;
t47 = t54 * pkin(7);
t31 = t54 * t43 + t47;
t50 = sin(qJ(4));
t53 = cos(qJ(4));
t12 = t53 * t30 - t50 * t31;
t60 = qJD(3) * t83;
t55 = cos(qJ(2));
t72 = pkin(1) * qJD(2);
t63 = t55 * t72;
t22 = t51 * t60 + t54 * t63;
t23 = -t51 * t63 + t54 * t60;
t2 = qJD(4) * t12 + t53 * t22 + t50 * t23;
t13 = t50 * t30 + t53 * t31;
t3 = -qJD(4) * t13 - t50 * t22 + t53 * t23;
t94 = t3 * mrSges(5,1) - t2 * mrSges(5,2);
t86 = -pkin(7) - pkin(6);
t41 = t86 * t51;
t42 = t54 * pkin(6) + t47;
t20 = t53 * t41 - t50 * t42;
t61 = qJD(3) * t86;
t35 = t51 * t61;
t36 = t54 * t61;
t10 = qJD(4) * t20 + t53 * t35 + t50 * t36;
t21 = t50 * t41 + t53 * t42;
t11 = -qJD(4) * t21 - t50 * t35 + t53 * t36;
t93 = t11 * mrSges(5,1) - t10 * mrSges(5,2);
t92 = t96 * t55;
t91 = qJD(3) + qJD(4);
t39 = -t54 * mrSges(4,1) + t51 * mrSges(4,2);
t90 = (t96 * mrSges(4,3) - mrSges(3,2)) * t55 + (t39 - mrSges(3,1)) * t52;
t89 = 2 * m(5);
t32 = -t50 * t51 + t53 * t54;
t17 = t91 * t32;
t33 = t50 * t54 + t53 * t51;
t18 = t91 * t33;
t6 = t18 * mrSges(5,1) + t17 * mrSges(5,2);
t88 = 0.2e1 * t6;
t19 = -t32 * mrSges(5,1) + t33 * mrSges(5,2);
t87 = 0.2e1 * t19;
t84 = t55 * pkin(1);
t80 = Ifges(4,6) * t51;
t73 = Ifges(5,5) * t17 - Ifges(5,6) * t18;
t71 = pkin(3) * qJD(4);
t70 = qJD(3) * t51;
t69 = qJD(3) * t54;
t68 = qJD(4) * t50;
t67 = qJD(4) * t53;
t66 = 2 * mrSges(5,3);
t65 = t53 * t17 * mrSges(5,3);
t64 = mrSges(5,3) * t71;
t62 = pkin(3) * t70;
t45 = -t54 * pkin(3) - pkin(2);
t58 = mrSges(4,1) * t51 + mrSges(4,2) * t54;
t57 = t53 * t32 * t64 + Ifges(4,5) * t69 + t73 + (-mrSges(5,3) * pkin(3) * t18 + t33 * t64) * t50;
t56 = 0.2e1 * t33 * t17 * Ifges(5,1) - 0.2e1 * t32 * Ifges(5,2) * t18 + (-0.2e1 * Ifges(4,4) * t51 + t95 * t54) * t70 + (0.2e1 * Ifges(4,4) * t54 + t95 * t51) * t69 + 0.2e1 * (t32 * t17 - t33 * t18) * Ifges(5,4);
t44 = -pkin(2) - t84;
t38 = t45 - t84;
t37 = t52 * t72 + t62;
t34 = t58 * qJD(3);
t29 = (-mrSges(5,1) * t50 - mrSges(5,2) * t53) * t71;
t1 = [(t12 * t3 + t13 * t2 + t38 * t37) * t89 + t37 * t87 + t38 * t88 + 0.2e1 * t44 * t34 + (-t12 * t17 - t13 * t18 + t2 * t32 - t3 * t33) * t66 + 0.2e1 * (m(4) * (t92 * t43 + t44 * t52) + t90) * t72 + t56; m(5) * (t10 * t13 + t11 * t12 + t21 * t2 + t20 * t3 + t45 * t37 + t38 * t62) + (t45 + t38) * t6 + (t44 - pkin(2)) * t34 + (t37 + t62) * t19 + (m(4) * (-pkin(2) * t52 + t92 * pkin(6)) + t90) * t72 + ((-t11 - t3) * t33 + (t10 + t2) * t32 - (t13 + t21) * t18 + (-t12 - t20) * t17) * mrSges(5,3) + t56; t62 * t87 + t45 * t88 - 0.2e1 * pkin(2) * t34 + (t21 * t10 + t20 * t11 + t45 * t62) * t89 + (t10 * t32 - t11 * t33 - t20 * t17 - t21 * t18) * t66 + t56; -t58 * t63 + (t39 * t43 - t80) * qJD(3) + (-t65 + m(5) * (-t12 * t68 + t13 * t67 + t2 * t50 + t3 * t53)) * pkin(3) + t57 + t94; (t39 * pkin(6) - t80) * qJD(3) + (-t65 + m(5) * (t10 * t50 + t11 * t53 - t20 * t68 + t21 * t67)) * pkin(3) + t57 + t93; 0.2e1 * t29; t73 + t94; t73 + t93; t29; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

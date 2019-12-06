% Calculate time derivative of joint inertia matrix for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:42
% EndTime: 2019-12-05 17:07:43
% DurationCPUTime: 0.59s
% Computational Cost: add. (829->128), mult. (1964->193), div. (0->0), fcn. (1469->6), ass. (0->71)
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t99 = t52 ^ 2 + t55 ^ 2;
t51 = sin(qJ(5));
t54 = cos(qJ(5));
t34 = t51 * t55 + t54 * t52;
t92 = qJD(4) + qJD(5);
t19 = t92 * t34;
t33 = -t51 * t52 + t54 * t55;
t98 = t33 * t19;
t18 = t92 * t33;
t97 = t34 * t18;
t96 = Ifges(5,1) - Ifges(5,2);
t53 = sin(qJ(3));
t44 = t53 * pkin(2) + pkin(7);
t84 = -pkin(8) - t44;
t31 = t84 * t52;
t48 = t55 * pkin(8);
t32 = t55 * t44 + t48;
t12 = t54 * t31 - t51 * t32;
t61 = qJD(4) * t84;
t56 = cos(qJ(3));
t73 = pkin(2) * qJD(3);
t64 = t56 * t73;
t23 = t52 * t61 + t55 * t64;
t24 = -t52 * t64 + t55 * t61;
t2 = t12 * qJD(5) + t54 * t23 + t51 * t24;
t13 = t51 * t31 + t54 * t32;
t3 = -t13 * qJD(5) - t51 * t23 + t54 * t24;
t95 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t87 = -pkin(8) - pkin(7);
t42 = t87 * t52;
t43 = t55 * pkin(7) + t48;
t21 = t54 * t42 - t51 * t43;
t62 = qJD(4) * t87;
t36 = t52 * t62;
t37 = t55 * t62;
t10 = t21 * qJD(5) + t54 * t36 + t51 * t37;
t22 = t51 * t42 + t54 * t43;
t11 = -t22 * qJD(5) - t51 * t36 + t54 * t37;
t94 = t11 * mrSges(6,1) - t10 * mrSges(6,2);
t93 = t99 * t56;
t59 = t52 * mrSges(5,1) + t55 * mrSges(5,2);
t35 = t59 * qJD(4);
t40 = -t55 * mrSges(5,1) + t52 * mrSges(5,2);
t91 = (t99 * mrSges(5,3) - mrSges(4,2)) * t56 + (t40 - mrSges(4,1)) * t53;
t90 = 2 * m(6);
t6 = t19 * mrSges(6,1) + t18 * mrSges(6,2);
t89 = 0.2e1 * t6;
t20 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t88 = 0.2e1 * t20;
t85 = t56 * pkin(2);
t81 = Ifges(5,6) * t52;
t74 = Ifges(6,5) * t18 - Ifges(6,6) * t19;
t72 = pkin(4) * qJD(5);
t71 = qJD(4) * t52;
t70 = qJD(4) * t55;
t69 = qJD(5) * t51;
t68 = qJD(5) * t54;
t67 = 2 * mrSges(6,3);
t66 = t54 * t18 * mrSges(6,3);
t65 = mrSges(6,3) * t72;
t63 = pkin(4) * t71;
t46 = -t55 * pkin(4) - pkin(3);
t58 = t54 * t33 * t65 + Ifges(5,5) * t70 + t74 + (-mrSges(6,3) * pkin(4) * t19 + t34 * t65) * t51;
t57 = 0.2e1 * Ifges(6,1) * t97 - 0.2e1 * Ifges(6,2) * t98 + (-0.2e1 * Ifges(5,4) * t52 + t96 * t55) * t71 + (0.2e1 * Ifges(5,4) * t55 + t96 * t52) * t70 + 0.2e1 * (t33 * t18 - t34 * t19) * Ifges(6,4);
t45 = -pkin(3) - t85;
t39 = t46 - t85;
t38 = t53 * t73 + t63;
t30 = (-mrSges(6,1) * t51 - mrSges(6,2) * t54) * t72;
t1 = [(t97 - t98) * t90; m(6) * (-t12 * t19 + t13 * t18 + t2 * t34 + t3 * t33); t38 * t88 + t39 * t89 + 0.2e1 * t45 * t35 + (t12 * t3 + t13 * t2 + t39 * t38) * t90 + (-t12 * t18 - t13 * t19 + t2 * t33 - t3 * t34) * t67 + 0.2e1 * (m(5) * (t93 * t44 + t45 * t53) + t91) * t73 + t57; m(6) * (t10 * t34 + t11 * t33 + t22 * t18 - t21 * t19); m(6) * (t10 * t13 + t11 * t12 + t22 * t2 + t21 * t3 + t46 * t38 + t39 * t63) + (t46 + t39) * t6 + (t45 - pkin(3)) * t35 + (t38 + t63) * t20 + (m(5) * (-pkin(3) * t53 + t93 * pkin(7)) + t91) * t73 + ((-t11 - t3) * t34 + (t10 + t2) * t33 - (t13 + t22) * t19 + (-t12 - t21) * t18) * mrSges(6,3) + t57; t63 * t88 + t46 * t89 - 0.2e1 * pkin(3) * t35 + (t22 * t10 + t21 * t11 + t46 * t63) * t90 + (t10 * t33 - t11 * t34 - t21 * t18 - t22 * t19) * t67 + t57; -t35 + m(6) * (t18 * t51 - t19 * t54 + (-t33 * t51 + t34 * t54) * qJD(5)) * pkin(4) - t6; -t59 * t64 + (t40 * t44 - t81) * qJD(4) + (-t66 + m(6) * (-t12 * t69 + t13 * t68 + t2 * t51 + t3 * t54)) * pkin(4) + t58 + t95; (t40 * pkin(7) - t81) * qJD(4) + (-t66 + m(6) * (t10 * t51 + t11 * t54 - t21 * t69 + t22 * t68)) * pkin(4) + t58 + t94; 0.2e1 * t30; -t6; t74 + t95; t74 + t94; t30; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

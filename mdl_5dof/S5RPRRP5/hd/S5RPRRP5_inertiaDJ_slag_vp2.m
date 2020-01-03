% Calculate time derivative of joint inertia matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (314->76), mult. (761->105), div. (0->0), fcn. (435->6), ass. (0->42)
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t73 = t38 ^ 2 + t40 ^ 2;
t72 = 2 * Ifges(5,4) - 2 * Ifges(6,5);
t28 = cos(pkin(8)) * pkin(1) + pkin(2);
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t65 = pkin(1) * sin(pkin(8));
t50 = t41 * t28 - t39 * t65;
t7 = t50 * qJD(3);
t71 = t73 * t7;
t70 = Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3);
t57 = t39 * t28 + t41 * t65;
t54 = qJD(4) * t40;
t69 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t73) * t7;
t46 = t38 * mrSges(6,1) - t40 * mrSges(6,3);
t18 = t46 * qJD(4);
t68 = 0.2e1 * t18;
t10 = pkin(7) + t57;
t67 = t71 * t10;
t66 = t71 * pkin(7);
t55 = qJD(4) * t38;
t53 = qJD(5) * t40;
t52 = m(6) * t53;
t23 = -t40 * mrSges(5,1) + t38 * mrSges(5,2);
t8 = t57 * qJD(3);
t51 = (-mrSges(4,1) + t23) * t8;
t49 = mrSges(6,2) * t53 + Ifges(6,6) * t55 + (Ifges(6,4) + Ifges(5,5)) * t54;
t47 = t38 * mrSges(5,1) + t40 * mrSges(5,2);
t22 = -t40 * mrSges(6,1) - t38 * mrSges(6,3);
t45 = -t40 * pkin(4) - t38 * qJ(5);
t21 = -pkin(3) + t45;
t11 = -pkin(4) * t55 + qJ(5) * t54 + t38 * qJD(5);
t44 = t45 * mrSges(6,2) - Ifges(5,6) * t38;
t43 = (-t38 * t72 + t70 * t40) * t55 + (t70 * t38 + t72 * t40) * t54;
t42 = m(6) * t45 + t22 + t23;
t30 = mrSges(6,2) * t54;
t19 = t47 * qJD(4);
t9 = -pkin(3) - t50;
t4 = t21 - t50;
t3 = t8 - t11;
t1 = [t43 + 0.2e1 * t51 + 0.2e1 * t69 + 0.2e1 * m(6) * (t4 * t3 + t67) + 0.2e1 * m(5) * (t9 * t8 + t67) + 0.2e1 * m(4) * (-t50 * t8 + t57 * t7) + 0.2e1 * t3 * t22 + t4 * t68 + 0.2e1 * t9 * t19; 0; 0; t51 + (t3 - t11) * t22 + (t9 - pkin(3)) * t19 + (t4 + t21) * t18 + m(6) * (-t11 * t4 + t21 * t3 + t66) + m(5) * (-pkin(3) * t8 + t66) + t69 + t43; 0; -0.2e1 * pkin(3) * t19 + t21 * t68 + 0.2e1 * (-m(6) * t21 - t22) * t11 + t43; t10 * t52 + (m(6) * (-pkin(4) * t38 + qJ(5) * t40) - t46 - t47) * t7 + (t42 * t10 + t44) * qJD(4) + t49; m(6) * t11 + ((-mrSges(5,2) + mrSges(6,3)) * t40 + (-mrSges(5,1) - mrSges(6,1)) * t38) * qJD(4); t44 * qJD(4) + (t42 * qJD(4) + t52) * pkin(7) + t49; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t30 + (t10 * t54 + t38 * t7) * m(6); m(6) * t55; m(6) * pkin(7) * t54 + t30; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

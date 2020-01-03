% Calculate joint inertia matrix for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR6_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:39
% DurationCPUTime: 0.34s
% Computational Cost: add. (341->105), mult. (685->155), div. (0->0), fcn. (580->6), ass. (0->44)
t59 = pkin(3) + qJ(2);
t39 = sin(pkin(7));
t41 = cos(pkin(7));
t58 = t39 ^ 2 + t41 ^ 2;
t42 = sin(qJ(5));
t43 = cos(qJ(5));
t38 = sin(pkin(8));
t54 = t38 * t41;
t14 = t39 * t43 + t42 * t54;
t15 = t39 * t42 - t43 * t54;
t57 = Ifges(6,5) * t15 + Ifges(6,6) * t14;
t32 = t38 ^ 2;
t40 = cos(pkin(8));
t34 = t40 ^ 2;
t56 = m(4) + m(5) * (t34 + t32);
t18 = t39 * mrSges(5,1) + mrSges(5,3) * t54;
t5 = -t14 * mrSges(6,1) + t15 * mrSges(6,2);
t55 = t18 - t5;
t53 = t40 * t41;
t49 = -t39 * qJ(3) - pkin(1);
t16 = (-pkin(2) - qJ(4)) * t41 + t49;
t22 = t59 * t39;
t7 = t40 * t16 + t38 * t22;
t52 = t58 * qJ(2) ^ 2;
t23 = t59 * t41;
t51 = t42 ^ 2 + t43 ^ 2;
t50 = Ifges(6,3) * t53 + t57;
t10 = (pkin(4) * t40 + pkin(6) * t38) * t41 + t23;
t4 = t39 * pkin(6) + t7;
t1 = t43 * t10 - t42 * t4;
t2 = t42 * t10 + t43 * t4;
t47 = -t1 * t42 + t2 * t43;
t46 = -mrSges(6,1) * t42 - mrSges(6,2) * t43;
t6 = -t38 * t16 + t40 * t22;
t19 = -t39 * mrSges(5,2) - mrSges(5,3) * t53;
t8 = -mrSges(6,2) * t53 + t14 * mrSges(6,3);
t9 = mrSges(6,1) * t53 - t15 * mrSges(6,3);
t45 = -t42 * t9 + t43 * t8 + t19;
t30 = t41 * mrSges(4,2);
t29 = t39 * mrSges(3,2);
t24 = mrSges(5,1) * t53;
t21 = -t41 * pkin(2) + t49;
t3 = -t39 * pkin(4) - t6;
t11 = [Ifges(6,1) * t15 ^ 2 - 0.2e1 * pkin(1) * t29 + 0.2e1 * t1 * t9 + 0.2e1 * t6 * t18 + 0.2e1 * t7 * t19 + 0.2e1 * t2 * t8 + 0.2e1 * t21 * t30 + 0.2e1 * t23 * t24 + 0.2e1 * t3 * t5 + Ifges(2,3) + (0.2e1 * Ifges(6,4) * t15 + Ifges(6,2) * t14) * t14 + (-0.2e1 * t21 * mrSges(4,3) + (Ifges(4,2) + Ifges(3,1) + Ifges(5,3)) * t39) * t39 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t23 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(3) * (pkin(1) ^ 2 + t52) + m(4) * (t21 ^ 2 + t52) + (0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * t23 * mrSges(5,2) * t38 + (t32 * Ifges(5,1) + Ifges(3,2) + Ifges(4,3)) * t41 + ((0.2e1 * Ifges(5,4) * t38 + Ifges(5,2) * t40) * t41 + t50 + t57) * t40 + 0.2e1 * (-Ifges(5,5) * t38 - Ifges(5,6) * t40 + Ifges(3,4) + Ifges(4,6)) * t39) * t41 + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * qJ(2) * t58; -m(3) * pkin(1) - t41 * mrSges(3,1) - t39 * mrSges(4,3) + t29 + t30 - t55 * t38 + t45 * t40 + m(6) * (t38 * t3 + t47 * t40) + m(5) * (-t38 * t6 + t40 * t7) + m(4) * t21; m(3) + m(6) * (t51 * t34 + t32) + t56; t55 * t40 + (m(4) * qJ(2) + mrSges(4,1)) * t39 + t45 * t38 + m(6) * (-t40 * t3 + t47 * t38) + m(5) * (t38 * t7 + t40 * t6); m(6) * (-0.1e1 + t51) * t40 * t38; m(6) * (t51 * t32 + t34) + t56; -mrSges(5,2) * t54 + t42 * t8 + t43 * t9 + t24 + m(6) * (t43 * t1 + t42 * t2) + m(5) * t23; 0; 0; m(6) * t51 + m(5); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t50; t46 * t40; t46 * t38; t43 * mrSges(6,1) - t42 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;

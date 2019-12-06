% Calculate joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:00
% EndTime: 2019-12-05 18:22:02
% DurationCPUTime: 0.33s
% Computational Cost: add. (313->104), mult. (563->131), div. (0->0), fcn. (373->6), ass. (0->44)
t38 = cos(qJ(4));
t61 = t38 ^ 2;
t36 = sin(qJ(4));
t50 = t36 ^ 2 + t61;
t60 = 0.2e1 * t50;
t27 = t36 * mrSges(6,2);
t16 = -t38 * mrSges(6,1) + t27;
t59 = 0.2e1 * t16;
t52 = t36 * mrSges(5,2);
t17 = -t38 * mrSges(5,1) + t52;
t58 = 0.2e1 * t17;
t57 = m(6) * pkin(4);
t37 = sin(qJ(2));
t56 = pkin(1) * t37;
t55 = t38 * pkin(4);
t39 = cos(qJ(2));
t25 = t39 * pkin(1) + pkin(2);
t34 = sin(pkin(8));
t35 = cos(pkin(8));
t7 = t35 * t25 - t34 * t56;
t54 = t7 * mrSges(4,1);
t8 = t34 * t25 + t35 * t56;
t53 = t8 * mrSges(4,2);
t51 = t36 * mrSges(6,3);
t49 = 0.2e1 * mrSges(6,3);
t24 = -t35 * pkin(2) - pkin(3);
t23 = t34 * pkin(2) + pkin(7);
t48 = t50 * t23;
t47 = (Ifges(5,6) + Ifges(6,6)) * t38 + (Ifges(5,5) + Ifges(6,5)) * t36;
t5 = -pkin(3) - t7;
t46 = t35 * mrSges(4,1) - t34 * mrSges(4,2);
t45 = -mrSges(5,1) * t36 - mrSges(5,2) * t38;
t44 = mrSges(5,3) * t60;
t43 = Ifges(3,3) + Ifges(4,3) + (Ifges(6,2) + Ifges(5,2)) * t61 + ((Ifges(6,1) + Ifges(5,1)) * t36 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t38) * t36;
t42 = (t39 * mrSges(3,1) - t37 * mrSges(3,2)) * pkin(1);
t26 = t38 * qJ(5);
t15 = t24 - t55;
t10 = t38 * t23 + t26;
t9 = (-qJ(5) - t23) * t36;
t6 = pkin(7) + t8;
t3 = t5 - t55;
t2 = t38 * t6 + t26;
t1 = (-qJ(5) - t6) * t36;
t4 = [t6 * t44 + m(5) * (t50 * t6 ^ 2 + t5 ^ 2) + (-t1 * t36 + t2 * t38) * t49 + m(3) * (t37 ^ 2 + t39 ^ 2) * pkin(1) ^ 2 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + t43 + Ifges(2,3) + 0.2e1 * t42 + 0.2e1 * t54 - 0.2e1 * t53 + t3 * t59 + t5 * t58; t54 - t53 + (t5 + t24) * t17 + (t15 + t3) * t16 + t42 + m(6) * (t9 * t1 + t10 * t2 + t15 * t3) + m(5) * (t24 * t5 + t6 * t48) + (m(4) * (t34 * t8 + t35 * t7) + t46) * pkin(2) + ((t10 + t2) * t38 + (-t1 - t9) * t36) * mrSges(6,3) + (t50 * t6 + t48) * mrSges(5,3) + t43; t15 * t59 + t24 * t58 + (t10 * t38 - t9 * t36) * t49 + t23 * t44 + m(6) * (t10 ^ 2 + t15 ^ 2 + t9 ^ 2) + m(5) * (t50 * t23 ^ 2 + t24 ^ 2) + t43 + (0.2e1 * t46 + m(4) * (t34 ^ 2 + t35 ^ 2) * pkin(2)) * pkin(2); m(6) * (t1 * t38 + t2 * t36); m(6) * (t10 * t36 + t9 * t38); m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t60; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t45 * t6 + (m(6) * t1 - t51) * pkin(4) + t47; t9 * mrSges(6,1) - t10 * mrSges(6,2) + t45 * t23 + (m(6) * t9 - t51) * pkin(4) + t47; -t52 - t27 + (mrSges(5,1) + mrSges(6,1) + t57) * t38; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t57) * pkin(4); m(6) * t3 + t16; m(6) * t15 + t16; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;

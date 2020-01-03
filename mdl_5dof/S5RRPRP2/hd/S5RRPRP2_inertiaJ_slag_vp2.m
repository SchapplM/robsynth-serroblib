% Calculate joint inertia matrix for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:27
% EndTime: 2019-12-31 19:49:28
% DurationCPUTime: 0.32s
% Computational Cost: add. (294->89), mult. (543->105), div. (0->0), fcn. (336->6), ass. (0->42)
t39 = cos(qJ(4));
t34 = t39 ^ 2;
t37 = sin(qJ(4));
t54 = t37 ^ 2 + t34;
t67 = -m(6) * pkin(4) - mrSges(6,1);
t66 = (mrSges(5,3) + mrSges(6,2)) * t54;
t15 = -t39 * mrSges(6,1) - t37 * mrSges(6,3);
t65 = 0.2e1 * t15;
t16 = -t39 * mrSges(5,1) + t37 * mrSges(5,2);
t64 = 0.2e1 * t16;
t35 = sin(pkin(8));
t25 = t35 * pkin(2) + pkin(7);
t40 = cos(qJ(2));
t27 = t40 * pkin(1) + pkin(2);
t36 = cos(pkin(8));
t38 = sin(qJ(2));
t60 = pkin(1) * t38;
t10 = t35 * t27 + t36 * t60;
t8 = pkin(7) + t10;
t63 = t54 * t25 * t8;
t62 = t54 * t8 ^ 2;
t61 = m(6) * t37;
t58 = t36 * pkin(2);
t9 = t36 * t27 - t35 * t60;
t57 = t9 * mrSges(4,1);
t56 = t10 * mrSges(4,2);
t29 = t37 * mrSges(6,2);
t55 = t54 * t25 ^ 2;
t53 = qJ(5) * t39;
t51 = t36 * mrSges(4,1) - t35 * mrSges(4,2);
t50 = t39 * pkin(4) + t37 * qJ(5);
t49 = -pkin(3) - t50;
t48 = (t40 * mrSges(3,1) - t38 * mrSges(3,2)) * pkin(1);
t47 = Ifges(3,3) + Ifges(4,3) + (Ifges(6,3) + Ifges(5,2)) * t34 + ((Ifges(6,1) + Ifges(5,1)) * t37 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t39) * t37;
t45 = mrSges(6,2) * t53 - pkin(4) * t29 + (Ifges(5,6) - Ifges(6,6)) * t39 + (Ifges(6,4) + Ifges(5,5)) * t37;
t44 = 0.2e1 * t66;
t43 = m(6) * t53 + (mrSges(6,3) - mrSges(5,2)) * t39 + (-mrSges(5,1) + t67) * t37;
t26 = -pkin(3) - t58;
t11 = t49 - t58;
t7 = -pkin(3) - t9;
t3 = t49 - t9;
t1 = [0.2e1 * t57 - 0.2e1 * t56 + t3 * t65 + t7 * t64 + Ifges(2,3) + 0.2e1 * t48 + t44 * t8 + m(6) * (t3 ^ 2 + t62) + m(5) * (t7 ^ 2 + t62) + m(4) * (t10 ^ 2 + t9 ^ 2) + m(3) * (t38 ^ 2 + t40 ^ 2) * pkin(1) ^ 2 + t47; t57 - t56 + (t7 + t26) * t16 + (t11 + t3) * t15 + t48 + m(6) * (t11 * t3 + t63) + m(5) * (t26 * t7 + t63) + (m(4) * (t10 * t35 + t36 * t9) + t51) * pkin(2) + t47 + (t25 + t8) * t66; t11 * t65 + t26 * t64 + m(6) * (t11 ^ 2 + t55) + m(5) * (t26 ^ 2 + t55) + t44 * t25 + t47 + (0.2e1 * t51 + m(4) * (t35 ^ 2 + t36 ^ 2) * pkin(2)) * pkin(2); 0; 0; m(4) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t54; t43 * t8 + t45; t43 * t25 + t45; m(6) * t50 - t15 - t16; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); t8 * t61 + t29; t25 * t61 + t29; -m(6) * t39; t67; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

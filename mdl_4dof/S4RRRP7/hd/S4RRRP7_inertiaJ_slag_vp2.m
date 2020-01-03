% Calculate joint inertia matrix for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:09
% EndTime: 2019-12-31 17:20:11
% DurationCPUTime: 0.39s
% Computational Cost: add. (255->133), mult. (539->180), div. (0->0), fcn. (375->4), ass. (0->53)
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t65 = t41 ^ 2 + t43 ^ 2;
t64 = 2 * pkin(5);
t44 = cos(qJ(2));
t63 = pkin(5) * t44;
t62 = Ifges(4,4) * t41;
t61 = Ifges(4,4) * t43;
t60 = Ifges(5,5) * t41;
t59 = Ifges(5,5) * t43;
t58 = Ifges(4,6) * t44;
t57 = Ifges(5,6) * t44;
t42 = sin(qJ(2));
t56 = t41 * t42;
t55 = t42 * t43;
t18 = -t44 * pkin(2) - t42 * pkin(6) - pkin(1);
t54 = t43 * t18;
t53 = Ifges(5,2) + Ifges(4,3);
t15 = t44 * mrSges(5,1) + mrSges(5,2) * t55;
t4 = t41 * t18 + t43 * t63;
t52 = t65 * pkin(6) ^ 2;
t51 = -Ifges(5,6) * t56 + (-Ifges(5,4) - Ifges(4,5)) * t55;
t49 = t41 * mrSges(4,1) + t43 * mrSges(4,2);
t48 = t41 * mrSges(5,1) - t43 * mrSges(5,3);
t47 = -pkin(3) * t41 + qJ(4) * t43;
t46 = pkin(5) ^ 2;
t40 = t44 ^ 2;
t38 = t42 ^ 2;
t35 = t38 * t46;
t33 = Ifges(5,4) * t41;
t32 = Ifges(4,5) * t41;
t31 = Ifges(4,6) * t43;
t24 = Ifges(4,1) * t41 + t61;
t23 = Ifges(5,1) * t41 - t59;
t22 = Ifges(4,2) * t43 + t62;
t21 = -Ifges(5,3) * t43 + t60;
t20 = -t43 * mrSges(4,1) + t41 * mrSges(4,2);
t19 = -t43 * mrSges(5,1) - t41 * mrSges(5,3);
t17 = -t43 * pkin(3) - t41 * qJ(4) - pkin(2);
t16 = -mrSges(5,2) * t56 - t44 * mrSges(5,3);
t14 = -t44 * mrSges(4,1) - mrSges(4,3) * t55;
t13 = t44 * mrSges(4,2) - mrSges(4,3) * t56;
t11 = t49 * t42;
t10 = t48 * t42;
t9 = (pkin(5) - t47) * t42;
t8 = -Ifges(4,5) * t44 + (Ifges(4,1) * t43 - t62) * t42;
t7 = -Ifges(5,4) * t44 + (Ifges(5,1) * t43 + t60) * t42;
t6 = -t58 + (-Ifges(4,2) * t41 + t61) * t42;
t5 = -t57 + (Ifges(5,3) * t41 + t59) * t42;
t3 = -t41 * t63 + t54;
t2 = -t54 + (pkin(5) * t41 + pkin(3)) * t44;
t1 = -t44 * qJ(4) + t4;
t12 = [0.2e1 * t1 * t16 + 0.2e1 * t9 * t10 + 0.2e1 * t4 * t13 + 0.2e1 * t3 * t14 + 0.2e1 * t2 * t15 + Ifges(2,3) + (t38 + t40) * mrSges(3,3) * t64 + m(5) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2 + t35) + m(3) * (pkin(1) ^ 2 + t40 * t46 + t35) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t53) * t44 + t51) * t44 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t42 + 0.2e1 * Ifges(3,4) * t44 + t11 * t64 + (t7 + t8) * t43 + (t5 - t6 + t58) * t41) * t42; -pkin(2) * t11 + t9 * t19 + (m(5) * t9 + t10) * t17 + (-pkin(5) * mrSges(3,2) + Ifges(3,6) - t32 / 0.2e1 - t31 / 0.2e1 - t33 / 0.2e1) * t44 + (t1 * mrSges(5,2) + t4 * mrSges(4,3) - t5 / 0.2e1 + t6 / 0.2e1 + t57 / 0.2e1) * t43 + (t2 * mrSges(5,2) - t3 * mrSges(4,3) + t7 / 0.2e1 + t8 / 0.2e1) * t41 + ((t13 + t16) * t43 + (-t14 + t15) * t41 + m(5) * (t1 * t43 + t2 * t41) + m(4) * (-t3 * t41 + t4 * t43)) * pkin(6) + (Ifges(3,5) + (t23 / 0.2e1 + t24 / 0.2e1) * t43 + (t21 / 0.2e1 - t22 / 0.2e1) * t41 + (-m(4) * pkin(2) - mrSges(3,1) + t20) * pkin(5)) * t42; -0.2e1 * pkin(2) * t20 + 0.2e1 * t17 * t19 + Ifges(3,3) + (-t21 + t22) * t43 + (t23 + t24) * t41 + m(5) * (t17 ^ 2 + t52) + m(4) * (pkin(2) ^ 2 + t52) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(6) * t65; -Ifges(4,6) * t56 - pkin(3) * t15 + m(5) * (-pkin(3) * t2 + qJ(4) * t1) + qJ(4) * t16 + t1 * mrSges(5,3) - t2 * mrSges(5,1) - t4 * mrSges(4,2) + t3 * mrSges(4,1) - t53 * t44 - t51; -Ifges(5,6) * t43 + t31 + t32 + t33 + t47 * mrSges(5,2) + (m(5) * t47 - t48 - t49) * pkin(6); 0.2e1 * pkin(3) * mrSges(5,1) + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t53; m(5) * t2 + t15; (m(5) * pkin(6) + mrSges(5,2)) * t41; -m(5) * pkin(3) - mrSges(5,1); m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t12(1), t12(2), t12(4), t12(7); t12(2), t12(3), t12(5), t12(8); t12(4), t12(5), t12(6), t12(9); t12(7), t12(8), t12(9), t12(10);];
Mq = res;

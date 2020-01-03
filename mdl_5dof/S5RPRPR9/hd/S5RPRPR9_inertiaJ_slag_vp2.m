% Calculate joint inertia matrix for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:23:44
% DurationCPUTime: 0.37s
% Computational Cost: add. (294->118), mult. (552->155), div. (0->0), fcn. (372->6), ass. (0->53)
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t46 = t31 ^ 2 + t33 ^ 2;
t44 = m(6) * t46;
t63 = m(5) + t44;
t32 = sin(qJ(3));
t25 = t32 ^ 2;
t34 = cos(qJ(3));
t27 = t34 ^ 2;
t62 = t25 + t27;
t61 = 0.2e1 * t62;
t60 = m(5) * pkin(3);
t59 = t46 * mrSges(6,3) - mrSges(5,2);
t49 = t34 * mrSges(6,3);
t10 = t32 * mrSges(6,1) + t31 * t49;
t35 = -pkin(3) - pkin(7);
t30 = cos(pkin(8));
t20 = -t30 * pkin(1) - pkin(2);
t21 = t32 * qJ(4);
t42 = t20 - t21;
t3 = t35 * t34 + t42;
t29 = sin(pkin(8));
t19 = t29 * pkin(1) + pkin(6);
t55 = pkin(4) + t19;
t8 = t55 * t32;
t1 = -t31 * t3 + t33 * t8;
t2 = t33 * t3 + t31 * t8;
t40 = t33 * t1 + t31 * t2;
t11 = -t32 * mrSges(6,2) - t33 * t49;
t52 = t31 * t11;
t58 = m(6) * t40 + t33 * t10 + t52;
t57 = -t31 / 0.2e1;
t56 = t34 * pkin(3);
t54 = Ifges(6,4) * t31;
t53 = Ifges(6,4) * t33;
t50 = t33 * t35;
t12 = t31 * mrSges(6,1) + t33 * mrSges(6,2);
t48 = t12 + mrSges(5,3);
t47 = t62 * t19 ^ 2;
t43 = t46 * t35;
t39 = t33 * mrSges(6,1) - t31 * mrSges(6,2);
t38 = -Ifges(6,5) * t31 - Ifges(6,6) * t33;
t36 = qJ(4) ^ 2;
t23 = Ifges(6,5) * t33;
t22 = Ifges(6,3) * t32;
t14 = Ifges(6,1) * t33 - t54;
t13 = -Ifges(6,2) * t31 + t53;
t9 = t55 * t34;
t7 = t39 * t34;
t6 = t42 - t56;
t5 = Ifges(6,5) * t32 + (-Ifges(6,1) * t31 - t53) * t34;
t4 = Ifges(6,6) * t32 + (-Ifges(6,2) * t33 - t54) * t34;
t15 = [0.2e1 * t1 * t10 + 0.2e1 * t2 * t11 + 0.2e1 * t9 * t7 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t20 * mrSges(4,2) - 0.2e1 * t6 * mrSges(5,3) + t22 + (Ifges(5,2) + Ifges(4,1)) * t32) * t32 + (0.2e1 * t6 * mrSges(5,2) - 0.2e1 * t20 * mrSges(4,1) - t31 * t5 - t33 * t4 + (Ifges(5,3) + Ifges(4,2)) * t34 + ((2 * Ifges(4,4)) + (2 * Ifges(5,6)) + t38) * t32) * t34 + m(6) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + m(5) * (t6 ^ 2 + t47) + m(4) * (t20 ^ 2 + t47) + (mrSges(5,1) + mrSges(4,3)) * t19 * t61 + (0.2e1 * t30 * mrSges(3,1) - 0.2e1 * t29 * mrSges(3,2) + m(3) * (t29 ^ 2 + t30 ^ 2) * pkin(1)) * pkin(1); (m(6) * t9 + t7) * t32 - t58 * t34; m(3) + m(6) * (t46 * t27 + t25) + (m(4) / 0.2e1 + m(5) / 0.2e1) * t61; qJ(4) * t7 + m(6) * (qJ(4) * t9 + t40 * t35) + t9 * t12 + t33 * t5 / 0.2e1 + t4 * t57 + t35 * t52 + t10 * t50 - t40 * mrSges(6,3) + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,6) * t57 + t23 / 0.2e1 - pkin(3) * mrSges(5,1)) * t32 + (-Ifges(5,5) + Ifges(4,6) + qJ(4) * mrSges(5,1) + t14 * t57 - t33 * t13 / 0.2e1) * t34 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t34 + (-mrSges(4,1) + mrSges(5,2) - t60) * t32) * t19; (mrSges(4,1) + t59) * t34 + m(5) * (t21 + t56) + m(6) * (-t34 * t43 + t21) + (-mrSges(4,2) + t48) * t32; -0.2e1 * pkin(3) * mrSges(5,2) - t31 * t13 + t33 * t14 + Ifges(5,1) + Ifges(4,3) + m(6) * (t46 * t35 ^ 2 + t36) + m(5) * (pkin(3) ^ 2 + t36) - 0.2e1 * mrSges(6,3) * t43 + 0.2e1 * t48 * qJ(4); (m(5) * t19 + mrSges(5,1)) * t32 + t58; -t63 * t34; t35 * t44 - t59 - t60; t63; t1 * mrSges(6,1) - t2 * mrSges(6,2) + t34 * t38 + t22; -t7; mrSges(6,1) * t50 + t23 + (-mrSges(6,2) * t35 - Ifges(6,6)) * t31; t39; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;

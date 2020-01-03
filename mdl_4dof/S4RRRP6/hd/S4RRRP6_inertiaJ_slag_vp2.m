% Calculate joint inertia matrix for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:11
% DurationCPUTime: 0.45s
% Computational Cost: add. (251->126), mult. (526->176), div. (0->0), fcn. (381->4), ass. (0->53)
t69 = Ifges(4,5) + Ifges(5,5);
t57 = Ifges(4,6) + Ifges(5,6);
t68 = 2 * pkin(5);
t67 = -2 * mrSges(5,3);
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t66 = t69 * t43 + t57 * t45;
t46 = cos(qJ(2));
t64 = pkin(5) * t46;
t63 = Ifges(4,4) * t43;
t62 = Ifges(4,4) * t45;
t61 = Ifges(5,4) * t43;
t60 = Ifges(5,4) * t45;
t44 = sin(qJ(2));
t59 = t43 * t44;
t58 = t44 * t45;
t56 = Ifges(4,3) + Ifges(5,3);
t55 = -qJ(4) - pkin(6);
t9 = mrSges(5,1) * t59 + mrSges(5,2) * t58;
t54 = t69 * t58;
t18 = -t46 * pkin(2) - t44 * pkin(6) - pkin(1);
t4 = t43 * t18 + t45 * t64;
t51 = t43 ^ 2 + t45 ^ 2;
t50 = qJ(4) * t44;
t20 = -t45 * mrSges(5,1) + t43 * mrSges(5,2);
t49 = mrSges(4,1) * t43 + mrSges(4,2) * t45;
t48 = pkin(5) ^ 2;
t42 = t46 ^ 2;
t40 = t44 ^ 2;
t38 = t40 * t48;
t32 = -t45 * pkin(3) - pkin(2);
t26 = Ifges(4,1) * t43 + t62;
t25 = Ifges(5,1) * t43 + t60;
t24 = Ifges(4,2) * t45 + t63;
t23 = Ifges(5,2) * t45 + t61;
t22 = t55 * t45;
t21 = -t45 * mrSges(4,1) + t43 * mrSges(4,2);
t19 = t55 * t43;
t17 = (pkin(3) * t43 + pkin(5)) * t44;
t16 = -t46 * mrSges(4,1) - mrSges(4,3) * t58;
t15 = -t46 * mrSges(5,1) - mrSges(5,3) * t58;
t14 = t46 * mrSges(4,2) - mrSges(4,3) * t59;
t13 = t46 * mrSges(5,2) - mrSges(5,3) * t59;
t12 = t45 * t18;
t10 = t49 * t44;
t8 = -Ifges(4,5) * t46 + (Ifges(4,1) * t45 - t63) * t44;
t7 = -Ifges(5,5) * t46 + (Ifges(5,1) * t45 - t61) * t44;
t6 = -Ifges(4,6) * t46 + (-Ifges(4,2) * t43 + t62) * t44;
t5 = -Ifges(5,6) * t46 + (-Ifges(5,2) * t43 + t60) * t44;
t3 = -t43 * t64 + t12;
t2 = -t43 * t50 + t4;
t1 = -t45 * t50 + t12 + (-pkin(5) * t43 - pkin(3)) * t46;
t11 = [0.2e1 * t1 * t15 + 0.2e1 * t2 * t13 + 0.2e1 * t4 * t14 + 0.2e1 * t3 * t16 + 0.2e1 * t17 * t9 + Ifges(2,3) + (t40 + t42) * mrSges(3,3) * t68 + m(5) * (t1 ^ 2 + t17 ^ 2 + t2 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2 + t38) + m(3) * (pkin(1) ^ 2 + t42 * t48 + t38) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t56) * t46 - t54) * t46 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t44 + 0.2e1 * Ifges(3,4) * t46 + t10 * t68 + (t7 + t8) * t45 + (t57 * t46 - t5 - t6) * t43) * t44; m(5) * (t19 * t1 + t32 * t17 - t22 * t2) - pkin(2) * t10 + t19 * t15 + t17 * t20 - t22 * t13 + t32 * t9 + (t2 * mrSges(5,3) + t4 * mrSges(4,3) + t5 / 0.2e1 + t6 / 0.2e1 + (m(4) * t4 + t14) * pkin(6)) * t45 + (-t1 * mrSges(5,3) - t3 * mrSges(4,3) + t7 / 0.2e1 + t8 / 0.2e1 + (-m(4) * t3 - t16) * pkin(6)) * t43 + (Ifges(3,5) + (-m(4) * pkin(2) - mrSges(3,1) + t21) * pkin(5) + (t25 / 0.2e1 + t26 / 0.2e1) * t45 + (-t23 / 0.2e1 - t24 / 0.2e1) * t43) * t44 + (Ifges(3,6) - pkin(5) * mrSges(3,2) - t66 / 0.2e1) * t46; -0.2e1 * pkin(2) * t21 + 0.2e1 * t32 * t20 + Ifges(3,3) + 0.2e1 * t51 * pkin(6) * mrSges(4,3) + m(5) * (t19 ^ 2 + t22 ^ 2 + t32 ^ 2) + m(4) * (t51 * pkin(6) ^ 2 + pkin(2) ^ 2) + (t22 * t67 + t23 + t24) * t45 + (t19 * t67 + t25 + t26) * t43; t3 * mrSges(4,1) + t1 * mrSges(5,1) - t4 * mrSges(4,2) - t2 * mrSges(5,2) - t56 * t46 - t57 * t59 + (m(5) * t1 + t15) * pkin(3) + t54; t19 * mrSges(5,1) + t22 * mrSges(5,2) - t49 * pkin(6) + (m(5) * t19 - t43 * mrSges(5,3)) * pkin(3) + t66; (m(5) * pkin(3) + 0.2e1 * mrSges(5,1)) * pkin(3) + t56; m(5) * t17 + t9; m(5) * t32 + t20; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1), t11(2), t11(4), t11(7); t11(2), t11(3), t11(5), t11(8); t11(4), t11(5), t11(6), t11(9); t11(7), t11(8), t11(9), t11(10);];
Mq = res;

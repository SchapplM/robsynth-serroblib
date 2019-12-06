% Calculate joint inertia matrix for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:38
% EndTime: 2019-12-05 16:41:39
% DurationCPUTime: 0.24s
% Computational Cost: add. (170->71), mult. (347->81), div. (0->0), fcn. (172->4), ass. (0->29)
t28 = cos(qJ(4));
t25 = t28 ^ 2;
t26 = sin(qJ(4));
t41 = t26 ^ 2 + t25;
t50 = -m(6) * pkin(4) - mrSges(6,1);
t49 = (mrSges(5,3) + mrSges(6,2)) * t41;
t8 = -t28 * mrSges(6,1) - t26 * mrSges(6,3);
t48 = 0.2e1 * t8;
t27 = sin(qJ(3));
t15 = t27 * pkin(2) + pkin(7);
t47 = t41 * pkin(7) * t15;
t46 = m(6) * t26;
t29 = cos(qJ(3));
t44 = t29 * pkin(2);
t18 = t26 * mrSges(6,2);
t43 = t41 * t15 ^ 2;
t42 = t41 * pkin(7) ^ 2;
t40 = qJ(5) * t28;
t9 = -t28 * mrSges(5,1) + t26 * mrSges(5,2);
t38 = t28 * pkin(4) + t26 * qJ(5);
t5 = -pkin(3) - t38;
t37 = Ifges(4,3) + (Ifges(6,3) + Ifges(5,2)) * t25 + ((Ifges(6,1) + Ifges(5,1)) * t26 + 0.2e1 * (Ifges(5,4) - Ifges(6,5)) * t28) * t26;
t36 = (t29 * mrSges(4,1) - t27 * mrSges(4,2)) * pkin(2);
t34 = mrSges(6,2) * t40 - pkin(4) * t18 + (Ifges(5,6) - Ifges(6,6)) * t28 + (Ifges(6,4) + Ifges(5,5)) * t26;
t33 = 0.2e1 * t49;
t32 = m(6) * t40 + (mrSges(6,3) - mrSges(5,2)) * t28 + (-mrSges(5,1) + t50) * t26;
t16 = -pkin(3) - t44;
t1 = t5 - t44;
t2 = [m(2) + m(3) + m(4) + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t41; 0; t1 * t48 + 0.2e1 * t16 * t9 + Ifges(3,3) + 0.2e1 * t36 + m(6) * (t1 ^ 2 + t43) + m(5) * (t16 ^ 2 + t43) + m(4) * (t27 ^ 2 + t29 ^ 2) * pkin(2) ^ 2 + t33 * t15 + t37; 0; (t16 - pkin(3)) * t9 + (t1 + t5) * t8 + t36 + m(6) * (t5 * t1 + t47) + m(5) * (-pkin(3) * t16 + t47) + t37 + (pkin(7) + t15) * t49; -0.2e1 * pkin(3) * t9 + t5 * t48 + m(6) * (t5 ^ 2 + t42) + m(5) * (pkin(3) ^ 2 + t42) + t33 * pkin(7) + t37; m(6) * t38 - t8 - t9; t15 * t32 + t34; pkin(7) * t32 + t34; Ifges(6,2) + Ifges(5,3) + 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * qJ(5) * mrSges(6,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2); -m(6) * t28; t15 * t46 + t18; pkin(7) * t46 + t18; t50; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;

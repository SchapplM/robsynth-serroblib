% Calculate joint inertia matrix for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:32
% EndTime: 2019-12-31 16:51:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (152->59), mult. (265->80), div. (0->0), fcn. (145->4), ass. (0->24)
t18 = cos(qJ(4));
t32 = t18 ^ 2;
t16 = sin(qJ(4));
t6 = -t18 * mrSges(5,1) + t16 * mrSges(5,2);
t31 = -mrSges(4,1) + t6;
t27 = t16 ^ 2 + t32;
t25 = t27 * mrSges(5,3);
t30 = -0.2e1 * t6;
t17 = sin(qJ(3));
t19 = cos(qJ(3));
t20 = -pkin(1) - pkin(2);
t4 = -t17 * qJ(2) + t19 * t20;
t29 = t4 * mrSges(4,1);
t5 = t19 * qJ(2) + t17 * t20;
t28 = t5 * mrSges(4,2);
t3 = -pkin(6) + t5;
t26 = t27 * t3;
t24 = t27 * t17;
t23 = -mrSges(5,1) * t16 - mrSges(5,2) * t18;
t22 = Ifges(5,2) * t32 + Ifges(4,3) + (Ifges(5,1) * t16 + 0.2e1 * Ifges(5,4) * t18) * t16;
t15 = t19 ^ 2;
t13 = t17 ^ 2;
t2 = pkin(3) - t4;
t1 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t29 + 0.2e1 * t28 + 0.2e1 * qJ(2) * mrSges(3,3) + t2 * t30 + Ifges(3,2) + Ifges(2,3) - 0.2e1 * t3 * t25 + m(5) * (t27 * t3 ^ 2 + t2 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + t22; -m(3) * pkin(1) - mrSges(3,1) + t31 * t19 + (mrSges(4,2) - t25) * t17 + m(5) * (-t19 * t2 + t3 * t24) + m(4) * (t17 * t5 + t19 * t4); m(3) + m(4) * (t13 + t15) + m(5) * (t27 * t13 + t15); m(5) * (-pkin(3) * t2 + pkin(6) * t26) + t29 - t28 + (pkin(3) + t2) * t6 + (-t27 * pkin(6) + t26) * mrSges(5,3) - t22; -t17 * mrSges(4,2) + (m(5) * pkin(6) + mrSges(5,3)) * t24 + (m(5) * pkin(3) - t31) * t19; pkin(3) * t30 + m(5) * (t27 * pkin(6) ^ 2 + pkin(3) ^ 2) + 0.2e1 * pkin(6) * t25 + t22; (-mrSges(5,2) * t3 - Ifges(5,6)) * t18 + (-mrSges(5,1) * t3 - Ifges(5,5)) * t16; t23 * t17; Ifges(5,5) * t16 + Ifges(5,6) * t18 + t23 * pkin(6); Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

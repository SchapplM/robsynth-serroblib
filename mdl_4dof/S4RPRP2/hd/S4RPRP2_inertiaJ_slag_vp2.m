% Calculate joint inertia matrix for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:23
% EndTime: 2018-11-14 13:49:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (81->39), mult. (110->42), div. (0->0), fcn. (51->2), ass. (0->16)
t12 = mrSges(4,2) + mrSges(5,2);
t10 = -pkin(1) - pkin(2);
t8 = sin(qJ(3));
t9 = cos(qJ(3));
t5 = t9 * qJ(2) + t8 * t10;
t17 = t12 * t5;
t16 = t12 * t8;
t15 = m(5) * pkin(3);
t4 = -t8 * qJ(2) + t9 * t10;
t14 = t4 * mrSges(4,1);
t13 = -mrSges(4,1) - mrSges(5,1);
t11 = Ifges(4,3) + Ifges(5,3);
t3 = t5 ^ 2;
t2 = -pkin(3) + t4;
t1 = t8 * t5;
t6 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t14 - 0.2e1 * t2 * mrSges(5,1) + 0.2e1 * qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + 0.2e1 * t17 + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + m(4) * (t4 ^ 2 + t3) + m(5) * (t2 ^ 2 + t3) + t11; -m(3) * pkin(1) - mrSges(3,1) + t13 * t9 + t16 + m(4) * (t9 * t4 + t1) + m(5) * (t9 * t2 + t1); m(3) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t8 ^ 2 + t9 ^ 2); t2 * t15 + t14 - t17 + (-pkin(3) + t2) * mrSges(5,1) - t11; -t16 + (-t13 + t15) * t9; (0.2e1 * mrSges(5,1) + t15) * pkin(3) + t11; 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1) t6(2) t6(4) t6(7); t6(2) t6(3) t6(5) t6(8); t6(4) t6(5) t6(6) t6(9); t6(7) t6(8) t6(9) t6(10);];
Mq  = res;

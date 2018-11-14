% Calculate joint inertia matrix for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:32
% EndTime: 2018-11-14 13:48:33
% DurationCPUTime: 0.08s
% Computational Cost: add. (72->36), mult. (124->42), div. (0->0), fcn. (64->4), ass. (0->15)
t17 = 2 * mrSges(5,3);
t8 = sin(pkin(6));
t16 = pkin(1) * t8;
t10 = sin(qJ(3));
t11 = cos(qJ(3));
t9 = cos(pkin(6));
t7 = pkin(1) * t9 + pkin(2);
t4 = t10 * t7 + t11 * t16;
t3 = -t10 * t16 + t11 * t7;
t15 = t3 * mrSges(4,1);
t14 = t4 * mrSges(4,2);
t13 = Ifges(5,2) + Ifges(4,3);
t2 = -pkin(3) - t3;
t1 = qJ(4) + t4;
t5 = [0.2e1 * t15 - 0.2e1 * t2 * mrSges(5,1) - 0.2e1 * t14 + t1 * t17 + Ifges(2,3) + Ifges(3,3) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t13 + (0.2e1 * t9 * mrSges(3,1) - 0.2e1 * t8 * mrSges(3,2) + m(3) * (t8 ^ 2 + t9 ^ 2) * pkin(1)) * pkin(1); 0; m(3) + m(4) + m(5); -t14 + t15 + m(5) * (-pkin(3) * t2 + qJ(4) * t1) + (t1 + qJ(4)) * mrSges(5,3) + (-t2 + pkin(3)) * mrSges(5,1) + t13; 0; 0.2e1 * pkin(3) * mrSges(5,1) + qJ(4) * t17 + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t13; m(5) * t2 - mrSges(5,1); 0; -m(5) * pkin(3) - mrSges(5,1); m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;

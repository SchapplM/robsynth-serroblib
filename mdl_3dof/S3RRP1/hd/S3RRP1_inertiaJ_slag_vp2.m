% Calculate joint inertia matrix for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3RRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_inertiaJ_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_inertiaJ_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_inertiaJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_inertiaJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRP1_inertiaJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:07
% DurationCPUTime: 0.08s
% Computational Cost: add. (35->25), mult. (64->31), div. (0->0), fcn. (15->2), ass. (0->8)
t8 = 2 * mrSges(4,3);
t7 = Ifges(4,2) + Ifges(3,3);
t3 = sin(qJ(2));
t4 = cos(qJ(2));
t6 = (t4 * mrSges(3,1) - t3 * mrSges(3,2)) * pkin(1);
t2 = -pkin(1) * t4 - pkin(2);
t1 = pkin(1) * t3 + qJ(3);
t5 = [-0.2e1 * t2 * mrSges(4,1) + t1 * t8 + Ifges(2,3) + 0.2e1 * t6 + m(4) * (t1 ^ 2 + t2 ^ 2) + m(3) * (t3 ^ 2 + t4 ^ 2) * pkin(1) ^ 2 + t7; m(4) * (-pkin(2) * t2 + qJ(3) * t1) + t6 + (t1 + qJ(3)) * mrSges(4,3) + (-t2 + pkin(2)) * mrSges(4,1) + t7; 0.2e1 * pkin(2) * mrSges(4,1) + qJ(3) * t8 + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t7; m(4) * t2 - mrSges(4,1); -m(4) * pkin(2) - mrSges(4,1); m(4);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t5(1) t5(2) t5(4); t5(2) t5(3) t5(5); t5(4) t5(5) t5(6);];
Mq  = res;

% Calculate joint inertia matrix for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_inertiaJ_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_inertiaJ_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_inertiaJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_inertiaJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPR1_inertiaJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:51
% EndTime: 2019-03-08 18:05:51
% DurationCPUTime: 0.07s
% Computational Cost: add. (36->22), mult. (55->27), div. (0->0), fcn. (24->2), ass. (0->9)
t3 = sin(qJ(3));
t4 = cos(qJ(3));
t5 = -pkin(1) - pkin(2);
t1 = -t3 * qJ(2) + t4 * t5;
t8 = t1 * mrSges(4,1);
t2 = t4 * qJ(2) + t3 * t5;
t7 = t2 * mrSges(4,2);
t6 = t4 * mrSges(4,1) - t3 * mrSges(4,2);
t9 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t8 + 0.2e1 * t7 + 0.2e1 * qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + m(4) * (t1 ^ 2 + t2 ^ 2); -mrSges(3,1) - m(3) * pkin(1) + m(4) * (t1 * t4 + t2 * t3) - t6; m(3) + m(4) * (t3 ^ 2 + t4 ^ 2); -Ifges(4,3) - t7 + t8; t6; Ifges(4,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t9(1) t9(2) t9(4); t9(2) t9(3) t9(5); t9(4) t9(5) t9(6);];
Mq  = res;

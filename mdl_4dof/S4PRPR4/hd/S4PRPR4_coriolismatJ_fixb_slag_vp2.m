% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (48->18), mult. (88->22), div. (0->0), fcn. (56->2), ass. (0->10)
t4 = sin(qJ(4));
t5 = cos(qJ(4));
t3 = t5 * mrSges(5,1) - t4 * mrSges(5,2);
t1 = -t5 ^ 2 * Ifges(5,4) + qJ(3) * t3 + (Ifges(5,4) * t4 + (-Ifges(5,1) + Ifges(5,2)) * t5) * t4;
t9 = t1 * qJD(2);
t7 = -t4 * mrSges(5,1) - t5 * mrSges(5,2);
t2 = mrSges(4,3) + (m(4) + m(5)) * qJ(3) - t7;
t8 = t2 * qJD(2);
t6 = -pkin(2) - pkin(5);
t10 = [0, 0, 0, -t3 * qJD(4); 0, t2 * qJD(3) + t1 * qJD(4), t8, t9 + ((-mrSges(5,2) * t6 - Ifges(5,6)) * t5 + (-mrSges(5,1) * t6 - Ifges(5,5)) * t4) * qJD(4); 0, -t8, 0, t7 * qJD(4); 0, -t9, 0, 0;];
Cq = t10;

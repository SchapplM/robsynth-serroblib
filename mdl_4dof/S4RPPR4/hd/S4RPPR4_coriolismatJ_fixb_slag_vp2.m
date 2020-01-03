% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (59->19), mult. (117->26), div. (0->0), fcn. (67->4), ass. (0->10)
t6 = sin(qJ(4));
t7 = cos(qJ(4));
t3 = t7 * mrSges(5,1) - t6 * mrSges(5,2);
t5 = sin(pkin(6)) * pkin(1) + qJ(3);
t1 = -t7 ^ 2 * Ifges(5,4) + t5 * t3 + (Ifges(5,4) * t6 + (-Ifges(5,1) + Ifges(5,2)) * t7) * t6;
t10 = t1 * qJD(1);
t8 = -t6 * mrSges(5,1) - t7 * mrSges(5,2);
t2 = mrSges(4,3) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t5 - t8;
t9 = t2 * qJD(1);
t4 = [t2 * qJD(3) + t1 * qJD(4), 0, t9, t10 + (-Ifges(5,5) * t6 - Ifges(5,6) * t7 + t8 * (-cos(pkin(6)) * pkin(1) - pkin(2) - pkin(5))) * qJD(4); 0, 0, 0, -t3 * qJD(4); -t9, 0, 0, t8 * qJD(4); -t10, 0, 0, 0;];
Cq = t4;

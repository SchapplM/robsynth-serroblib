% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:49
% EndTime: 2019-03-08 18:18:49
% DurationCPUTime: 0.07s
% Computational Cost: add. (62->16), mult. (158->26), div. (0->0), fcn. (145->4), ass. (0->11)
t11 = cos(pkin(5));
t14 = cos(qJ(2));
t8 = sin(pkin(5));
t9 = sin(qJ(2));
t4 = t11 * t9 + t8 * t14;
t15 = m(5) * t4;
t6 = t8 * pkin(2) + qJ(4);
t5 = m(5) * t6 + mrSges(5,3);
t12 = t5 * qJD(2);
t3 = -t11 * t14 + t8 * t9;
t1 = [0 (-t14 * mrSges(3,2) - t9 * mrSges(3,1) + t3 * mrSges(4,2) - t4 * mrSges(4,1) + m(4) * (-t11 * t4 - t3 * t8) * pkin(2) - t4 * mrSges(5,1) - t3 * mrSges(5,3) + m(5) * (-t6 * t3 + (-t11 * pkin(2) - pkin(3)) * t4)) * qJD(2) + qJD(4) * t15, 0, qJD(2) * t15; 0, t5 * qJD(4), 0, t12; 0, 0, 0, 0; 0, -t12, 0, 0;];
Cq  = t1;

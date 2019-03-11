% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR1
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
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:26
% EndTime: 2019-03-08 18:27:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (202->17), mult. (277->21), div. (0->0), fcn. (189->4), ass. (0->15)
t12 = sin(qJ(4));
t13 = cos(qJ(4));
t14 = t12 * mrSges(5,1) + t13 * mrSges(5,2);
t20 = t14 * qJD(4);
t10 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(3);
t11 = sin(pkin(6)) * pkin(1) + qJ(3);
t8 = t13 * t10 - t12 * t11;
t9 = t12 * t10 + t13 * t11;
t1 = t9 * mrSges(5,1) + t8 * mrSges(5,2);
t19 = t1 * qJD(1);
t18 = t1 * qJD(4);
t2 = mrSges(4,3) + m(4) * t11 + m(5) * (-t8 * t12 + t9 * t13) + t14;
t17 = t2 * qJD(1);
t16 = t14 * qJD(1);
t3 = [t2 * qJD(3) + t18, 0, t17, -t18 + t19; 0, 0, 0, 0; -t17 + t20, 0, 0, t16 - t20; -qJD(3) * t14 - t19, 0, -t16, 0;];
Cq  = t3;

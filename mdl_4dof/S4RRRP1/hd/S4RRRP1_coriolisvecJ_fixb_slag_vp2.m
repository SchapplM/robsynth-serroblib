% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:53
% EndTime: 2019-03-08 18:35:53
% DurationCPUTime: 0.26s
% Computational Cost: add. (333->55), mult. (851->81), div. (0->0), fcn. (396->4), ass. (0->35)
t41 = mrSges(4,1) + mrSges(5,1);
t23 = qJD(1) + qJD(2);
t27 = cos(qJ(2));
t39 = pkin(1) * qJD(1);
t19 = t23 * pkin(2) + t27 * t39;
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t25 = sin(qJ(2));
t34 = t25 * t39;
t15 = t24 * t19 + t26 * t34;
t43 = t24 * t25;
t31 = t26 * t27 - t43;
t17 = t31 * t39;
t37 = qJD(3) * t26;
t38 = qJD(3) * t24;
t29 = (t31 * qJD(2) - t25 * t38) * pkin(1);
t7 = qJD(1) * t29 + t19 * t37;
t49 = -t15 * t17 + (t15 * t37 + t24 * t7) * pkin(2);
t48 = mrSges(3,1) * t25 + mrSges(3,2) * t27;
t21 = t27 * pkin(1) + pkin(2);
t10 = t21 * t37 + t29;
t42 = t25 * t26;
t47 = t7 * (pkin(1) * t42 + t24 * t21) + t15 * t10;
t40 = mrSges(4,2) + mrSges(5,2);
t33 = -pkin(1) * t43 + t26 * t21;
t32 = -t24 * t27 - t42;
t28 = (t32 * qJD(2) - t25 * t37) * pkin(1);
t8 = qJD(1) * t28 - t19 * t38;
t30 = -t40 * t7 + t41 * t8;
t14 = t26 * t19 - t24 * t34;
t22 = qJD(3) + t23;
t16 = t32 * t39;
t12 = t22 * pkin(3) + t14;
t11 = -t21 * t38 + t28;
t1 = [m(4) * (t14 * t11 + t8 * t33 + t47) + m(5) * (t8 * (pkin(3) + t33) + t12 * t11 + t47) + (-t40 * t10 + t41 * t11) * t22 + t30 + t48 * pkin(1) * qJD(2) * (-qJD(1) - t23); (t40 * t17 - t41 * t16 + (-t41 * t24 - t40 * t26) * qJD(3) * pkin(2)) * t22 + t30 + t48 * t39 * (-qJD(2) + t23) + (t8 * (t26 * pkin(2) + pkin(3)) + (-pkin(2) * t38 - t16) * t12 + t49) * m(5) + ((-t14 * t38 + t8 * t26) * pkin(2) - t14 * t16 + t49) * m(4); (t40 * t14 + t41 * t15) * t22 + t30 + (t8 * pkin(3) + (t12 - t14) * t15) * m(5); 0;];
tauc  = t1(:);

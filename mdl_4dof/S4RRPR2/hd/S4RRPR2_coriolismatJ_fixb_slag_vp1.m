% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:31
% EndTime: 2019-07-18 18:16:33
% DurationCPUTime: 0.86s
% Computational Cost: add. (6549->82), mult. (5524->114), div. (0->0), fcn. (5452->6), ass. (0->71)
t174 = m(5) / 0.2e1;
t187 = m(4) / 0.2e1;
t124 = qJ(1) + qJ(2);
t122 = sin(t124);
t123 = cos(t124);
t149 = sin(qJ(1)) * pkin(1);
t119 = t123 * qJ(3);
t173 = pkin(2) + pkin(3);
t150 = sin(qJ(4));
t151 = cos(qJ(4));
t110 = -t122 * t150 - t123 * t151;
t111 = -t122 * t151 + t123 * t150;
t180 = -t111 * rSges(5,1) + t110 * rSges(5,2);
t184 = -t122 * t173 + t119 - t180;
t188 = t184 - t149;
t148 = cos(qJ(1)) * pkin(1);
t89 = t110 * rSges(5,1) + t111 * rSges(5,2);
t83 = t122 * qJ(3) + t123 * t173 - t89;
t78 = t83 + t148;
t44 = t78 * t122 + t188 * t123;
t45 = t83 * t122 + t123 * t184;
t152 = rSges(4,1) + pkin(2);
t104 = t123 * rSges(4,3) - t122 * t152 + t119;
t100 = t104 - t149;
t105 = t152 * t123 + (rSges(4,3) + qJ(3)) * t122;
t101 = t105 + t148;
t66 = t100 * t123 + t101 * t122;
t72 = t104 * t123 + t105 * t122;
t142 = (t45 + t44) * t174 + (t72 + t66) * t187;
t143 = (-t45 + t44) * t174 + (t66 - t72) * t187;
t4 = t143 - t142;
t195 = t4 * qJD(1);
t191 = -t78 * t180 + t188 * t89;
t194 = m(5) * qJD(1) * t191;
t182 = t89 * t122 + t123 * t180;
t193 = m(5) * t182;
t60 = -t193 / 0.2e1;
t61 = t193 / 0.2e1;
t136 = m(5) * qJD(4);
t34 = t83 * t180 - t89 * t184;
t192 = t191 * t136;
t29 = t184 * t78 - t188 * t83;
t168 = m(4) * (-t105 * t100 + t101 * t104);
t172 = m(3) * (t148 * (-t122 * rSges(3,1) - t123 * rSges(3,2)) + (t123 * rSges(3,1) - t122 * rSges(3,2)) * t149);
t178 = 0.4e1 * qJD(1);
t166 = m(4) * t66;
t165 = m(4) * t72;
t164 = m(5) * (t34 - t191);
t163 = m(5) * (-t34 - t191);
t159 = m(5) * t29;
t156 = m(5) * t44;
t155 = m(5) * t45;
t22 = 0.2e1 * t61;
t144 = t22 * qJD(3);
t137 = m(5) * qJD(2);
t130 = qJD(1) + qJD(2);
t33 = t155 + t165;
t31 = t156 + t166;
t21 = t60 + t61;
t20 = 0.2e1 * t60;
t18 = t21 * qJD(3);
t17 = t21 * qJD(4);
t16 = t20 * qJD(4);
t11 = t163 / 0.2e1;
t10 = t164 / 0.2e1;
t7 = t159 + t168 + t172;
t5 = t142 + t143;
t3 = t11 - t164 / 0.2e1;
t2 = t10 - t163 / 0.2e1;
t1 = t10 + t11;
t6 = [t7 * qJD(2) + t31 * qJD(3) - t192, t7 * qJD(1) + t5 * qJD(3) + t1 * qJD(4) + 0.2e1 * (t172 / 0.2e1 + t168 / 0.2e1 + t29 * t174) * qJD(2), t31 * qJD(1) + t5 * qJD(2) + t17, t1 * qJD(2) + t18 + t192 - t194; -t4 * qJD(3) + t2 * qJD(4) + (-t172 / 0.4e1 - t168 / 0.4e1 - t159 / 0.4e1) * t178, t33 * qJD(3) + t136 * t34, t33 * qJD(2) + t17 - t195, t2 * qJD(1) + t18 + (-t136 + t137) * t34; t4 * qJD(2) + t16 + (-t166 / 0.4e1 - t156 / 0.4e1) * t178, t195 + t16 + 0.4e1 * (-t165 / 0.4e1 - t155 / 0.4e1) * qJD(2), 0, t130 * t20 + t136 * t182; t3 * qJD(2) + t144 + t194, t3 * qJD(1) - t34 * t137 + t144, t130 * t22, 0;];
Cq  = t6;

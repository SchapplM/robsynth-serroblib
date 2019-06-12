% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-31 10:31
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-31 10:31:13
% EndTime: 2019-05-31 10:31:14
% DurationCPUTime: 0.28s
% Computational Cost: add. (264->54), mult. (390->92), div. (0->0), fcn. (295->8), ass. (0->53)
t103 = qJ(2) + qJ(3);
t100 = sin(t103);
t107 = cos(qJ(4));
t154 = r_i_i_C(1) * t107 + pkin(2);
t157 = t100 * t154;
t104 = sin(qJ(4));
t137 = qJD(4) * t107;
t101 = cos(t103);
t102 = qJD(2) + qJD(3);
t145 = t101 * t102;
t156 = t100 * t137 + t104 * t145;
t151 = pkin(5) + r_i_i_C(3);
t132 = t151 * t101;
t105 = sin(qJ(2));
t146 = pkin(1) * qJD(2);
t134 = t105 * t146;
t149 = pkin(2) * t100;
t155 = -t134 + (t132 - t149) * t102;
t138 = qJD(4) * t104;
t125 = t100 * t138;
t152 = r_i_i_C(1) * t125 + t156 * r_i_i_C(2);
t150 = pkin(1) * t105;
t147 = r_i_i_C(2) * t104;
t106 = sin(qJ(1));
t144 = t102 * t106;
t143 = t102 * t107;
t109 = cos(qJ(1));
t142 = t102 * t109;
t141 = t107 * t109;
t140 = qJD(1) * t106;
t139 = qJD(1) * t109;
t136 = t100 * t147;
t135 = qJD(1) * t147;
t133 = t151 * t100;
t131 = t151 * t106;
t130 = t100 * t143;
t120 = qJD(4) * t101 - qJD(1);
t119 = qJD(1) * t101 - qJD(4);
t118 = t152 * t109 + t140 * t157;
t117 = t154 * t101;
t116 = t154 * t109;
t115 = t109 * t100 * t135 + t152 * t106 + t139 * t132;
t114 = t120 * t104;
t113 = t100 * t142 + t119 * t106;
t108 = cos(qJ(2));
t112 = qJD(1) * (-pkin(1) * t108 - pkin(2) * t101 - t133);
t111 = -t108 * t146 + (-t117 - t133) * t102;
t110 = -t101 * r_i_i_C(2) * t137 + (-t101 * t138 - t130) * r_i_i_C(1) + t151 * t145 + (-t149 + t136) * t102;
t85 = -t119 * t141 + (t114 + t130) * t106;
t84 = t120 * t107 * t106 + (-t100 * t144 + t119 * t109) * t104;
t83 = t113 * t107 + t109 * t114;
t82 = t113 * t104 - t120 * t141;
t1 = [t85 * r_i_i_C(1) + t84 * r_i_i_C(2) - t155 * t106 + t109 * t112, (-t132 - t136 + t150) * t140 + t111 * t109 + t118, (-t106 * t135 - t151 * t142) * t100 + (-qJD(1) * t131 - t102 * t116) * t101 + t118, t82 * r_i_i_C(1) + t83 * r_i_i_C(2), 0; -t83 * r_i_i_C(1) + t82 * r_i_i_C(2) + t106 * t112 + t155 * t109, (-t150 - t157) * t139 + t111 * t106 + t115, -t117 * t144 + (-qJD(1) * t116 - t102 * t131) * t100 + t115, -t84 * r_i_i_C(1) + t85 * r_i_i_C(2), 0; 0, t110 - t134, t110, (-t101 * t143 + t125) * r_i_i_C(2) - t156 * r_i_i_C(1), 0;];
JaD_transl  = t1;

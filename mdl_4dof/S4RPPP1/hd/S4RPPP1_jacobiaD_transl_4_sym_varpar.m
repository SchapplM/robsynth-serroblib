% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPPP1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:58
% EndTime: 2019-02-26 19:29:58
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->24), mult. (126->36), div. (0->0), fcn. (112->6), ass. (0->25)
t141 = sin(pkin(4));
t160 = t141 * (pkin(3) + r_i_i_C(1) + qJ(2));
t159 = -r_i_i_C(2) - qJ(3);
t140 = sin(pkin(6));
t144 = sin(qJ(1));
t158 = t144 * t140;
t142 = cos(pkin(6));
t157 = t144 * t142;
t145 = cos(qJ(1));
t156 = t145 * t140;
t155 = t145 * t142;
t154 = qJD(1) * t144;
t153 = qJD(1) * t145;
t152 = t141 * qJD(2);
t151 = -pkin(2) - r_i_i_C(3) - qJ(4);
t149 = t140 * t154;
t148 = t142 * t153;
t143 = cos(pkin(4));
t147 = t143 * t157 + t156;
t146 = t143 * t156 + t157;
t137 = -t143 * t149 + t148;
t136 = t147 * qJD(1);
t135 = t146 * qJD(1);
t134 = -t143 * t148 + t149;
t1 = [-t146 * qJD(4) - (-t143 * t155 + t158) * qJD(3) + t145 * t152 + t159 * t136 + t151 * t137 + (-t145 * pkin(1) - t144 * t160) * qJD(1), t141 * t153, -t134, -t135; -(t143 * t158 - t155) * qJD(4) + t147 * qJD(3) + t144 * t152 + t159 * t134 + t151 * t135 + (-t144 * pkin(1) + t145 * t160) * qJD(1), t141 * t154, t136, t137; 0, 0, 0, 0;];
JaD_transl  = t1;

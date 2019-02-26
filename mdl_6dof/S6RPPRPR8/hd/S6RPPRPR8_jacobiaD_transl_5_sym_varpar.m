% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:45
% EndTime: 2019-02-26 20:29:45
% DurationCPUTime: 0.09s
% Computational Cost: add. (101->26), mult. (152->41), div. (0->0), fcn. (106->5), ass. (0->19)
t145 = pkin(9) + qJ(4);
t142 = sin(t145);
t143 = cos(t145);
t159 = r_i_i_C(3) + qJ(5);
t160 = pkin(4) - r_i_i_C(2);
t164 = -(t159 * t142 + t160 * t143) * qJD(4) + t143 * qJD(5);
t163 = t160 * qJD(4) - qJD(5);
t161 = t160 * t142 - t159 * t143 + pkin(3) * sin(pkin(9)) + qJ(2);
t148 = sin(qJ(1));
t158 = qJD(1) * t148;
t149 = cos(qJ(1));
t144 = qJD(1) * t149;
t157 = qJD(4) * t148;
t156 = qJD(4) * t149;
t154 = -pkin(1) - r_i_i_C(1) - pkin(7) - qJ(3);
t153 = qJD(1) * t160;
t151 = qJD(1) * t159;
t150 = qJD(2) - t164;
t1 = [-t148 * qJD(3) + t150 * t149 + (-t161 * t148 + t154 * t149) * qJD(1), t144, -t158 (t149 * t153 + t159 * t157) * t143 + (-t163 * t148 + t149 * t151) * t142, t142 * t157 - t143 * t144, 0; t149 * qJD(3) + t150 * t148 + (t154 * t148 + t161 * t149) * qJD(1), t158, t144 (t148 * t153 - t159 * t156) * t143 + (t148 * t151 + t163 * t149) * t142, -t142 * t156 - t143 * t158, 0; 0, 0, 0, t164, qJD(4) * t143, 0;];
JaD_transl  = t1;

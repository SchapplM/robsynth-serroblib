% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:23
% EndTime: 2019-02-26 19:45:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (59->23), mult. (202->41), div. (0->0), fcn. (192->8), ass. (0->23)
t162 = sin(pkin(11));
t165 = cos(pkin(11));
t168 = sin(qJ(2));
t169 = cos(qJ(2));
t161 = t168 * t162 - t165 * t169;
t159 = t161 * qJD(2);
t175 = pkin(3) - r_i_i_C(2);
t174 = r_i_i_C(3) + qJ(4);
t173 = qJD(2) * pkin(2);
t167 = cos(pkin(6));
t172 = t167 * t168;
t171 = t162 * t169 + t165 * t168;
t158 = t171 * t167;
t160 = t171 * qJD(2);
t170 = qJD(2) * t158;
t166 = cos(pkin(10));
t164 = sin(pkin(6));
t163 = sin(pkin(10));
t157 = t167 * t159;
t155 = t164 * t160;
t152 = t159 * t166 + t163 * t170;
t150 = t163 * t159 - t166 * t170;
t1 = [0 -(t158 * t163 + t161 * t166) * qJD(4) - t174 * (-t157 * t163 + t160 * t166) + t175 * t152 + (t163 * t172 - t166 * t169) * t173, 0, -t152, 0, 0; 0 -(-t158 * t166 + t161 * t163) * qJD(4) - t174 * (t157 * t166 + t160 * t163) + t175 * t150 + (-t163 * t169 - t166 * t172) * t173, 0, -t150, 0, 0; 0, -t175 * t155 + (t171 * qJD(4) - t174 * t159 - t168 * t173) * t164, 0, t155, 0, 0;];
JaD_transl  = t1;

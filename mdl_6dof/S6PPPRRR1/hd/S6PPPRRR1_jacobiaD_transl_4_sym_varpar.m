% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPPRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:54
% EndTime: 2019-02-26 19:38:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (65->32), mult. (208->71), div. (0->0), fcn. (252->14), ass. (0->40)
t181 = sin(pkin(12));
t190 = cos(pkin(6));
t204 = t181 * t190;
t182 = sin(pkin(8));
t191 = sin(qJ(4));
t203 = t182 * t191;
t192 = cos(qJ(4));
t202 = t182 * t192;
t183 = sin(pkin(7));
t184 = sin(pkin(6));
t201 = t183 * t184;
t200 = t183 * t190;
t189 = cos(pkin(7));
t199 = t184 * t189;
t186 = cos(pkin(13));
t198 = t186 * t189;
t187 = cos(pkin(12));
t197 = t187 * t190;
t188 = cos(pkin(8));
t196 = t188 * t191;
t195 = t188 * t192;
t180 = sin(pkin(13));
t175 = -t180 * t181 + t186 * t197;
t194 = t175 * t189 - t187 * t201;
t177 = -t180 * t187 - t186 * t204;
t193 = t177 * t189 + t181 * t201;
t185 = cos(pkin(14));
t179 = sin(pkin(14));
t178 = -t180 * t204 + t186 * t187;
t176 = t180 * t197 + t181 * t186;
t174 = -t186 * t201 + t189 * t190;
t173 = -t177 * t183 + t181 * t199;
t172 = -t175 * t183 - t187 * t199;
t171 = t184 * t180 * t185 + (t184 * t198 + t200) * t179;
t170 = t185 * t200 + (-t179 * t180 + t185 * t198) * t184;
t169 = t178 * t185 + t179 * t193;
t168 = -t178 * t179 + t185 * t193;
t167 = t176 * t185 + t179 * t194;
t166 = -t176 * t179 + t185 * t194;
t1 = [0, 0, 0 ((-t168 * t196 - t169 * t192 - t173 * t203) * r_i_i_C(1) + (-t168 * t195 + t169 * t191 - t173 * t202) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, 0 ((-t166 * t196 - t167 * t192 - t172 * t203) * r_i_i_C(1) + (-t166 * t195 + t167 * t191 - t172 * t202) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, 0 ((-t170 * t196 - t171 * t192 - t174 * t203) * r_i_i_C(1) + (-t170 * t195 + t171 * t191 - t174 * t202) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;

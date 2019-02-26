% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR9_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:13
% EndTime: 2019-02-26 21:33:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (118->36), mult. (337->52), div. (0->0), fcn. (312->6), ass. (0->27)
t175 = sin(pkin(6));
t193 = t175 * (pkin(3) + pkin(8) + r_i_i_C(1));
t192 = r_i_i_C(2) + qJ(3);
t177 = sin(qJ(2));
t178 = sin(qJ(1));
t191 = t178 * t177;
t179 = cos(qJ(2));
t190 = t178 * t179;
t180 = cos(qJ(1));
t189 = t180 * t177;
t188 = t180 * t179;
t187 = qJD(2) * t177;
t186 = qJD(2) * t179;
t184 = pkin(2) + r_i_i_C(3) + qJ(4);
t176 = cos(pkin(6));
t183 = t176 * t191;
t182 = t176 * t188;
t181 = qJD(2) * t176 + qJD(1);
t169 = t176 * t190 + t189;
t168 = t176 * t189 + t190;
t170 = t183 - t188;
t167 = -t182 + t191;
t166 = -qJD(1) * t183 - t178 * t187 + t181 * t188;
t165 = t169 * qJD(1) + t168 * qJD(2);
t164 = t168 * qJD(1) + t169 * qJD(2);
t163 = -qJD(1) * t182 - t180 * t186 + t181 * t191;
t1 = [-t167 * qJD(3) - t168 * qJD(4) - t192 * t165 - t184 * t166 + (-pkin(1) * t180 - t178 * t193) * qJD(1), -t170 * qJD(3) - t169 * qJD(4) + t184 * t163 - t192 * t164, -t163, -t164, 0, 0; t169 * qJD(3) - t170 * qJD(4) - t192 * t163 - t184 * t164 + (-pkin(1) * t178 + t180 * t193) * qJD(1), t168 * qJD(3) - t167 * qJD(4) - t184 * t165 + t192 * t166, t165, t166, 0, 0; 0 (qJD(3) * t177 + qJD(4) * t179 + (-t184 * t177 + t192 * t179) * qJD(2)) * t175, t175 * t187, t175 * t186, 0, 0;];
JaD_transl  = t1;

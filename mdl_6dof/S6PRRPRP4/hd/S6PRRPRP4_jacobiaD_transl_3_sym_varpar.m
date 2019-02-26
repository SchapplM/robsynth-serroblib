% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRP4_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP4_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobiaD_transl_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:00
% EndTime: 2019-02-26 20:03:00
% DurationCPUTime: 0.09s
% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
t179 = sin(pkin(10));
t181 = cos(pkin(10));
t184 = sin(qJ(2));
t182 = cos(pkin(6));
t186 = cos(qJ(2));
t194 = t182 * t186;
t200 = -t179 * t184 + t181 * t194;
t199 = pkin(8) + r_i_i_C(3);
t180 = sin(pkin(6));
t183 = sin(qJ(3));
t197 = t180 * t183;
t185 = cos(qJ(3));
t196 = t180 * t185;
t195 = t182 * t184;
t192 = t183 * r_i_i_C(1) + t185 * r_i_i_C(2);
t191 = t185 * r_i_i_C(1) - t183 * r_i_i_C(2) + pkin(2);
t175 = t179 * t186 + t181 * t195;
t190 = t179 * t194 + t181 * t184;
t189 = t179 * t195 - t181 * t186;
t188 = qJD(3) * t192;
t187 = qJD(2) * t191;
t172 = t190 * qJD(2);
t170 = t200 * qJD(2);
t1 = [0, -t199 * t172 + t189 * t187 + t190 * t188, t192 * t172 + ((-t179 * t197 + t185 * t189) * r_i_i_C(1) + (-t179 * t196 - t183 * t189) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t199 * t170 - t175 * t187 - t200 * t188, -t192 * t170 + ((-t175 * t185 + t181 * t197) * r_i_i_C(1) + (t175 * t183 + t181 * t196) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0 (-t186 * t188 + (-t191 * t184 + t199 * t186) * qJD(2)) * t180, -t192 * t186 * t180 * qJD(2) + ((-t182 * t183 - t184 * t196) * r_i_i_C(1) + (-t182 * t185 + t184 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
JaD_transl  = t1;

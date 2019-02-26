% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPRR7
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
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR7_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_transl_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:57
% EndTime: 2019-02-26 21:31:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
t179 = sin(pkin(6));
t198 = t179 * (pkin(8) + r_i_i_C(2));
t197 = pkin(2) + r_i_i_C(1);
t195 = r_i_i_C(3) + qJ(3);
t181 = sin(qJ(2));
t182 = sin(qJ(1));
t194 = t182 * t181;
t183 = cos(qJ(2));
t193 = t182 * t183;
t184 = cos(qJ(1));
t192 = t184 * t181;
t191 = t184 * t183;
t190 = qJD(2) * t181;
t180 = cos(pkin(6));
t189 = t180 * t194;
t188 = t180 * t191;
t187 = qJD(2) * t180 + qJD(1);
t186 = t180 * t193 + t192;
t185 = t180 * t192 + t193;
t174 = -qJD(1) * t189 - t182 * t190 + t187 * t191;
t173 = t186 * qJD(1) + t185 * qJD(2);
t172 = t185 * qJD(1) + t186 * qJD(2);
t171 = -qJD(1) * t188 - qJD(2) * t191 + t187 * t194;
t1 = [-(-t188 + t194) * qJD(3) - t197 * t174 - t195 * t173 + (-t184 * pkin(1) - t182 * t198) * qJD(1) -(t189 - t191) * qJD(3) - t195 * t172 + t197 * t171, -t171, 0, 0, 0; t186 * qJD(3) - t197 * t172 - t195 * t171 + (-t182 * pkin(1) + t184 * t198) * qJD(1), t185 * qJD(3) - t197 * t173 + t195 * t174, t173, 0, 0, 0; 0 (t181 * qJD(3) + (-t197 * t181 + t195 * t183) * qJD(2)) * t179, t179 * t190, 0, 0, 0;];
JaD_transl  = t1;

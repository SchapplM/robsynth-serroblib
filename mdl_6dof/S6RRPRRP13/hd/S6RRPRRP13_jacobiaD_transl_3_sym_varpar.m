% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP13_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP13_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobiaD_transl_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:54
% EndTime: 2019-02-26 21:52:54
% DurationCPUTime: 0.11s
% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
t175 = sin(pkin(6));
t194 = t175 * (pkin(8) + r_i_i_C(1));
t193 = pkin(2) - r_i_i_C(2);
t191 = r_i_i_C(3) + qJ(3);
t177 = sin(qJ(2));
t178 = sin(qJ(1));
t190 = t178 * t177;
t179 = cos(qJ(2));
t189 = t178 * t179;
t180 = cos(qJ(1));
t188 = t180 * t177;
t187 = t180 * t179;
t186 = qJD(2) * t177;
t176 = cos(pkin(6));
t185 = t176 * t190;
t184 = t176 * t187;
t183 = qJD(2) * t176 + qJD(1);
t182 = t176 * t189 + t188;
t181 = t176 * t188 + t189;
t170 = -qJD(1) * t185 - t178 * t186 + t183 * t187;
t169 = qJD(1) * t182 + qJD(2) * t181;
t168 = qJD(1) * t181 + qJD(2) * t182;
t167 = -qJD(1) * t184 - qJD(2) * t187 + t183 * t190;
t1 = [-(-t184 + t190) * qJD(3) - t193 * t170 - t191 * t169 + (-t180 * pkin(1) - t178 * t194) * qJD(1) -(t185 - t187) * qJD(3) - t191 * t168 + t193 * t167, -t167, 0, 0, 0; t182 * qJD(3) - t193 * t168 - t191 * t167 + (-t178 * pkin(1) + t180 * t194) * qJD(1), qJD(3) * t181 - t169 * t193 + t170 * t191, t169, 0, 0, 0; 0 (t177 * qJD(3) + (-t177 * t193 + t179 * t191) * qJD(2)) * t175, t175 * t186, 0, 0, 0;];
JaD_transl  = t1;

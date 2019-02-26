% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR7_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_transl_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:30
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.07s
% Computational Cost: add. (41->25), mult. (151->54), div. (0->0), fcn. (138->10), ass. (0->24)
t193 = r_i_i_C(3) + qJ(3);
t175 = sin(pkin(14));
t181 = cos(pkin(7));
t192 = t175 * t181;
t177 = sin(pkin(7));
t183 = sin(qJ(2));
t191 = t177 * t183;
t179 = cos(pkin(14));
t190 = t179 * t181;
t184 = cos(qJ(2));
t189 = t181 * t184;
t182 = cos(pkin(6));
t188 = t182 * t183;
t187 = t182 * t184;
t176 = sin(pkin(13));
t180 = cos(pkin(13));
t186 = -t176 * t184 - t180 * t188;
t185 = t176 * t188 - t180 * t184;
t178 = sin(pkin(6));
t174 = t185 * qJD(2);
t173 = (t176 * t187 + t180 * t183) * qJD(2);
t172 = t186 * qJD(2);
t171 = (t176 * t183 - t180 * t187) * qJD(2);
t1 = [0 (t173 * t192 + t174 * t179) * r_i_i_C(1) + (t173 * t190 - t174 * t175) * r_i_i_C(2) + t174 * pkin(2) + (-t185 * qJD(3) - t193 * t173) * t177, -t174 * t177, 0, 0, 0; 0 (t171 * t192 + t172 * t179) * r_i_i_C(1) + (t171 * t190 - t172 * t175) * r_i_i_C(2) + t172 * pkin(2) + (-t186 * qJD(3) - t193 * t171) * t177, -t172 * t177, 0, 0, 0; 0 (qJD(3) * t191 + ((-t175 * t189 - t179 * t183) * r_i_i_C(1) + (t175 * t183 - t179 * t189) * r_i_i_C(2) - t183 * pkin(2) + t193 * t184 * t177) * qJD(2)) * t178, t178 * qJD(2) * t191, 0, 0, 0;];
JaD_transl  = t1;
